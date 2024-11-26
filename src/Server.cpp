/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.30
*/
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>
#include <thread>
#include <cctype>
#include <cstdlib>
#include <random>
#include <set>
#include <utility>
#include <exception>
#include <signal.h>
#include <unistd.h>

#include <boost/asio.hpp>  
#include <boost/asio/thread_pool.hpp>
#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include <grpc/grpc.h>
#include <grpcpp/grpcpp.h>

#include "FedVectorDB.grpc.pb.h"

#include "database/VectorDB.hpp"
#include "index/BaseIndex.hpp"
#include "utils/DataType.hpp"
#include "utils/MetricType.hpp"
#include "utils/BenchLogger.hpp"
#include "utils/File_IO.h"
#include "middleware/BaseSiloConnector.hpp"
#include "middleware/FileSiloConnector.hpp"

using grpc::Server;
using grpc::ServerBuilder;
using grpc::ServerContext;
using grpc::ServerReader;
using grpc::ServerReaderWriter;
using grpc::ServerWriter;
using grpc::Status;
using grpc::Channel;
using grpc::ClientContext;
using grpc::ClientReader;
using grpc::ClientReaderWriter;
using grpc::ClientWriter;
using google::protobuf::Empty;

using FedVectorDB::KnnSql;
using FedVectorDB::CandSql;
using FedVectorDB::VectorDataCand;
using FedVectorDB::FedSqlService;

using Distance = EuclideanSquareDistance;


class SiloReceiver {
public:
    SiloReceiver(std::shared_ptr<grpc::Channel> channel, const int silo_id, const std::string& silo_ipaddr) 
        : m_stub_(FedSqlService::NewStub(channel)), m_silo_id(silo_id), m_silo_ipaddr(silo_ipaddr) {
        
        m_logger.Init();
    }

    void RequestKnnCand(const VectorDataType& query_data, const size_t& query_k) {
        ClientContext context;
        Empty response;
        KnnSql query;

        query.set_vid(query_data.vid);
        query.set_k(query_k);
        for (auto d : query_data.data) {
            query.add_data(d);
        }

        // std::cout << "Silo " << m_silo_id << ": ``RequestKnnCand``" << std::endl;

        Status status = m_stub_->ExecuteKnnQuery(&context, query, &response); 
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Request Knn candidates from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }

        float grpc_comm = query.ByteSizeLong() + response.ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);
    }

    void RequestPartialKnn(const size_t cand_num) {
        ClientContext context;
        VectorDataCand cand;
        CandSql query;

        query.set_k(cand_num);
        float grpc_comm = query.ByteSizeLong();

        std::unique_ptr<ClientReader<VectorDataCand> > reader(
            m_stub_->GetKnnCandidate(&context, query));

        m_partial_knn.clear();
        while (reader->Read(&cand)) {
            grpc_comm += cand.ByteSizeLong();

            const size_t dim = cand.data_size();
            VectorDataType knn_data(dim, cand.vid());
            for (int i=0; i<dim; ++i) {
                knn_data[i] = cand.data(i);
            }
            m_partial_knn.emplace_back(knn_data);
        }

        Status status = reader->Finish();
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Request partial Knn answers from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }

        m_logger.LogAddComm(grpc_comm);  
    }

    void ClearKnnCand() {
        ClientContext context;
        Empty request;
        Empty response;

        Status status = m_stub_->ClearKnnCandidate(&context, request, &response); 
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Clear Knn candidates from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }

        float grpc_comm = request.ByteSizeLong() + response.ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);
    } 

    std::vector<VectorDataType> GetPartialKnn() const {
        return m_partial_knn;
    }

    double GetQueryComm() const {
        return m_logger.GetQueryComm();
    }

    void InitBenchLogger() {
        m_logger.Init();
    }

    static void RequestKnnCandThread(SiloReceiver* silo_receiver, const VectorDataType& query_data, const size_t& query_k) {  
        silo_receiver->RequestKnnCand(query_data, query_k);  
    } 

    static void RequestPartialKnnThread(SiloReceiver* silo_receiver, const size_t cand_num) {  
        silo_receiver->RequestPartialKnn(cand_num);
    }  

    static void ClearKnnCandThread(SiloReceiver* silo_receiver) {  
        silo_receiver->ClearKnnCand();  
    }   

private:
    std::unique_ptr<FedSqlService::Stub> m_stub_;
    std::vector<VectorDataType> m_partial_knn;
    std::string m_silo_ipaddr;
    int m_silo_id;  
    BenchLogger m_logger;  
};

class FedSqlServer {
public:
    FedSqlServer(const std::string& silo_ip_filename) 
        : m_thread_num(std::thread::hardware_concurrency()), m_silo_ip_filename(silo_ip_filename) {

        ReadSiloIPaddr(silo_ip_filename, m_silo_ipaddr_list);
        if (m_silo_ipaddr_list.empty()) {
            throw std::invalid_argument("No data silo's IP addresses");  
            std::exit(EXIT_FAILURE);
        }

        std::cout << "Server is running\n";

        m_CreateSiloReceiver();
    }

    FedSqlServer(const std::string& silo_ip_filename, const std::string& query_filename, const std::string& output_filename="", const std::string& drf_option="equality", const int &query_k=128)
                : FedSqlServer(silo_ip_filename) {
        
        ProcessFedKnnQuery(query_filename, output_filename, drf_option, query_k);
    }

    FedSqlServer(const std::string& silo_ip_filename, const std::string& query_filename, const std::string& output_filename, const std::string& truth_filename, const std::string& drf_option, const int &query_k)
                : FedSqlServer(silo_ip_filename, query_filename, output_filename, drf_option, query_k) {
        
        if (!output_filename.empty() && !truth_filename.empty()) {
            std::cout << "Server's output file name: " << output_filename << std::endl;
            std::cout << "Server's ground truth file name: " << truth_filename << std::endl;
            EvaluateAnswer(output_filename, truth_filename);
        }
    }

    void ProcessFedKnnQuery(const std::string& query_filename, const std::string& output_filename="", const std::string& drf_option="equality", const int &query_k=128) {
        m_query_filename = query_filename;

        // Step 0: Initialize local variables
        m_InitBenchLogger();
        m_CheckDrfOption(drf_option);

        // Step 1: Read vector query from file
        ReadVectorQuery(query_filename, m_query_data_list, m_query_k_list);
        std::fill(m_query_k_list.begin(), m_query_k_list.end(), query_k);
        std::cout << "Server receives " << m_query_data_list.size() << " queries" << std::endl;

        // Step 2: Enumerate all vector queries
        size_t query_num = m_query_data_list.size();
        if (output_filename.empty()) {
            std::cout << "Server will dump query answer into std::cout" << std::endl;
            // Dump the result into std.cout by default
            for (size_t i=0; i<query_num; ++i) {
                std::vector<VectorDataType> answer = ProcessOneFedKnnQuery(m_query_data_list[i], m_query_k_list[i], drf_option);
            }
        } else {
            std::cout << "Server will dump query answer into file: " << output_filename << std::endl;
            // Dump the result into specific file
            std::vector<std::vector<VidType>> answer_list;
            for (size_t i=0; i<query_num; ++i) {
                std::vector<VectorDataType> answer = ProcessOneFedKnnQuery(m_query_data_list[i], m_query_k_list[i], drf_option);
                std::vector<VidType> answer_vid(answer.size());
                for (size_t j=0,sz=answer_vid.size(); j<sz; ++j) {
                    answer_vid[j] = answer[j].vid;
                }
                answer_list.emplace_back(answer_vid);
            }
            DumpGroundTruth(output_filename, answer_list);
            
            // std::ofstream file(output_filename);
            // if (!file.is_open()) {
            //     std::cerr << "Failed to open file for dumping vector data: " << output_filename << std::endl;
            //     std::exit(EXIT_FAILURE);
            // }
            // size_t dim  = (query_num==0) ? 1 : m_query_data_list[0].Dimension();
            // file << query_num << " " << dim << "\n";
            // for (size_t i=0; i<query_num; ++i) {
            //     std::vector<VectorDataType> answer = ProcessOneFedKnnQuery(m_query_data_list[i], m_query_k_list[i], drf_option);
            //     size_t answer_num = answer.size();
            //     file << answer_num << "\n";
            //     for (size_t j=0; j<answer_num; ++j) {
            //         file << answer[j].vid;
            //         for (size_t k=0; k<dim; ++k) {
            //             file << " " << answer[j].data[k];
            //         }
            //         file << "\n";
            //     }
            // }
            // file << m_logger.to_string();
            // file.close();
        }

        // Step 3: Print the log information
        m_logger.Print();
    }

    std::vector<VectorDataType> ProcessOneFedKnnQuery(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option="equality") {
        std::vector<VectorDataType> answer;
        
        if (m_IsBaselineMode(drf_option)) {
            answer = m_ProcessFedKnnQuery(query_data, query_k, drf_option);
        } else {
            if (drf_option == std::string("greedy"))
                answer = m_ProcessFedKnnQuery_ByGreedy(query_data, query_k, drf_option);
            else if (drf_option == std::string("greedyX"))
                answer = m_ProcessFedKnnQuery_ByGreedyX(query_data, query_k, drf_option);
        }

        return answer;
    }

    std::string to_string() const {
        std::stringstream ss;

        ss << "-------------- Service Log --------------\n";
        ss << m_logger.to_string();

        return ss.str();
    }

private:
    std::vector<VectorDataType> m_ProcessFedKnnQuery(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option) {
        m_logger.SetStartTimer();

        #ifdef LOCAL_DEBUG
        std::cout << "k = " << query_k << ", query_data = " << query_data.to_string() << std::endl;
        #endif

        // Step 0: Intialize some variables
        std::vector<VectorDataType> answer;

        // Step 1: Perform local knn query at each silo
        m_PerformLocalKnn(query_data, query_k);

        // Step 2: Determine the resource allocation for each silo
        std::vector<size_t> allocation = m_DetermineResourceAllocation(query_data, query_k, drf_option);
        #ifdef LOCAL_DEBUG
        std::cout << "Server's allocation: ";
        std::cout << "[";
        for (size_t i=0; i<allocation.size(); ++i) {
            std::cout << allocation[i];
            if (i+1 < allocation.size())
                std::cout <<", ";
        }
        std::cout << "]" << std::endl;
        #endif

        // Step 3: Retrieve specific amount of partial answers from each silo based on the allocation
        m_RetrievePartialKnn(allocation);
        if (m_IsNotPlaintext(drf_option)) { 
            for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
                if (allocation[silo_id] > 0) {
                    std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
                    answer.insert(answer.end(), partial_knn_list.begin(), partial_knn_list.end());
                    #ifdef LOCAL_DEBUG
                    assert(partial_knn_list.size() == allocation[silo_id]);
                    #endif
                }
            }
        } else {
            answer = m_PickGlobalKnn(query_data, query_k);
        }

        // Step 4: Clear Knn candidates at each silo
        m_ClearKnnCand();

        // Step 5: Print the log information
        m_logger.SetEndTimer();
        double query_comm = 0.0;
        for (int i=0; i<m_silo_num; ++i) {
            query_comm += m_silo_receiver_list[i]->GetQueryComm();
        }
        query_comm -= m_logger.GetQueryComm();
        double query_time = m_logger.GetDurationTime();
        m_logger.LogOneQuery(query_comm);

        std::cout << std::fixed << std::setprecision(6) 
                    << "Query #(" << query_data.vid << "): runtime = " << query_time/1000.0 << " [s], communication = " << query_comm/1024.0 << " [KB]" << std::endl;
    
        return answer;
    }

    std::vector<VectorDataType> m_ProcessFedKnnQuery_ByGreedyX(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option) {
        m_logger.SetStartTimer();
        Distance dist_function;

        // Step 0: Intialize some variables
        std::vector<VectorDataType> answer;

        // Step 1: Perform local knn query at each silo
        m_PerformLocalKnn(query_data, query_k);

        // Step 2: Determine the resource allocation for each silo
        std::vector<VectorDimensionType> dist_list(m_silo_num);
        std::set<std::pair<VectorDimensionType,int>> minQ;
        std::vector<size_t> allocation(m_silo_num);
        for (size_t rid=0, total_answer=0; total_answer<query_k; ++rid) {
            
            if (rid == 0) {
                // Step 3a: Retrieve one candidate from each silo
                std::fill(allocation.begin(), allocation.end(), (size_t)1);
                m_RetrievePartialKnn(allocation);
                for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
                    std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
                    answer.insert(answer.end(), partial_knn_list.begin(), partial_knn_list.end());
                    total_answer += allocation[silo_id];
                    VectorDataType& vector_data = partial_knn_list[0];
                    dist_list[silo_id] = dist_function(vector_data, query_data);
                    minQ.insert(std::make_pair(dist_list[silo_id], silo_id));
                }

            } else {

                // Step 3b: Retrieve one candidate from the most relevant silo
                int silo_id = minQ.begin()->second;
                allocation[silo_id] = 1;
                SiloReceiver::RequestPartialKnnThread(m_silo_receiver_list[silo_id].get(), allocation[silo_id]);
                std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
                answer.insert(answer.end(), partial_knn_list.begin(), partial_knn_list.end());
                total_answer += allocation[silo_id];
                VectorDataType& vector_data = partial_knn_list[0];
                minQ.erase(std::make_pair(dist_list[silo_id], silo_id));
                dist_list[silo_id] = dist_function(vector_data, query_data);
                minQ.insert(std::make_pair(dist_list[silo_id], silo_id));
            }

        }
        

        // Step 4: Clear Knn candidates at each silo
        m_ClearKnnCand();

        // Step 5: Print the log information
        m_logger.SetEndTimer();
        double query_comm = 0.0;
        for (int i=0; i<m_silo_num; ++i) {
            query_comm += m_silo_receiver_list[i]->GetQueryComm();
        }
        query_comm -= m_logger.GetQueryComm();
        double query_time = m_logger.GetDurationTime();
        m_logger.LogOneQuery(query_comm);

        std::cout << std::fixed << std::setprecision(6) 
                    << "Query #(" << query_data.vid << "): runtime = " << query_time/1000.0 << " [s], communication = " << query_comm/1024.0 << " [KB]" << std::endl;
    
        return answer;
    }

    std::vector<VectorDataType> m_ProcessFedKnnQuery_ByGreedy(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option) {
        m_logger.SetStartTimer();
        Distance dist_function;

        #ifdef LOCAL_DEBUG
        std::cout << "k = " << query_k << ", query_data = " << query_data.to_string() << std::endl;
        #endif

        // Step 0: Intialize some variables
        std::vector<VectorDataType> answer;

        // Step 1: Perform local knn query at each silo
        m_PerformLocalKnn(query_data, query_k);

        // Step 2: Determine the resource allocation for each silo
        std::vector<VectorDimensionType> dist_list(m_silo_num, -1e7);
        std::set<std::pair<VectorDimensionType,int>> minQ;
        for (size_t rid=0, total_answer=0; total_answer<query_k; ++rid) {
            std::vector<size_t> allocation = m_DetermineResourceAllocationByGreedy(query_data, query_k, drf_option, total_answer, rid, minQ);
            #ifdef LOCAL_DEBUG
            std::cout << "Server's allocation at round " << rid << ": ";
            std::cout << "[";
            for (size_t i=0; i<allocation.size(); ++i) {
                std::cout << allocation[i];
                if (i+1 < allocation.size())
                    std::cout <<", ";
            }
            std::cout << "]" << std::endl;
            #endif

            // Step 3: Retrieve specific amount of partial answers from each silo based on the allocation
            m_RetrievePartialKnn(allocation);
            for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
                if (allocation[silo_id] > 0) {
                    std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
                    answer.insert(answer.end(), partial_knn_list.begin(), partial_knn_list.end());
                    #ifdef LOCAL_DEBUG
                    assert(partial_knn_list.size() == allocation[silo_id]);
                    #endif
                    total_answer += allocation[silo_id];

                    VectorDimensionType max_dist = -1e7;
                    for (const auto& vector_data : partial_knn_list) {
                        VectorDimensionType dist = dist_function(vector_data, query_data);
                        max_dist = std::max(max_dist, dist);
                    }
                    if (max_dist > dist_list[silo_id]) {
                        minQ.erase(std::make_pair(dist_list[silo_id], silo_id));
                        dist_list[silo_id] = max_dist;
                        minQ.insert(std::make_pair(max_dist, silo_id));
                    }
                }
            }
            #ifdef LOCAL_DEBUG
            std::cout << "Server's minQ at round " << rid << ": ";
            std::cout << "minQ.size() = " << minQ.size() << ", ";
            std::cout << "minQ.top() = Silo #" << minQ.begin()->second << std::endl;
            #endif
        }
        

        // Step 4: Clear Knn candidates at each silo
        m_ClearKnnCand();

        // Step 5: Print the log information
        m_logger.SetEndTimer();
        double query_comm = 0.0;
        for (int i=0; i<m_silo_num; ++i) {
            query_comm += m_silo_receiver_list[i]->GetQueryComm();
        }
        query_comm -= m_logger.GetQueryComm();
        double query_time = m_logger.GetDurationTime();
        m_logger.LogOneQuery(query_comm);

        std::cout << std::fixed << std::setprecision(6) 
                    << "Query #(" << query_data.vid << "): runtime = " << query_time/1000.0 << " [s], communication = " << query_comm/1024.0 << " [KB]" << std::endl;
    
        return answer;
    }

    void m_ClearKnnCand() {
        boost::asio::thread_pool pool(m_thread_num);

        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            boost::asio::post(pool, std::bind(SiloReceiver::ClearKnnCandThread, m_silo_receiver_list[silo_id].get())); 
        }

        pool.join();
    }

    void m_RetrievePartialKnn(const std::vector<size_t>& allocation) {
        if (allocation.size() != m_silo_num) {
            std::cerr << "The scheme for dominant resource allocation only covers " << allocation.size() 
                        << " out of " << m_silo_num << " data silos" << std::endl;
            std::exit(EXIT_FAILURE);   
        }

        boost::asio::thread_pool pool(m_thread_num);

        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            if (allocation[silo_id] > 0) {
                boost::asio::post(pool, std::bind(SiloReceiver::RequestPartialKnnThread, m_silo_receiver_list[silo_id].get(), allocation[silo_id])); 
            }
        }

        pool.join();
    }

    void m_PerformLocalKnn(const VectorDataType& query_data, const size_t& query_k) {
        boost::asio::thread_pool pool(m_thread_num);

        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            boost::asio::post(pool, std::bind(SiloReceiver::RequestKnnCandThread, m_silo_receiver_list[silo_id].get(), query_data, query_k)); 
        }

        pool.join();
    }

    std::vector<size_t> m_DetermineResourceAllocation(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option) {
        std::string input_drf_option(drf_option);
        std::transform(input_drf_option.begin(), input_drf_option.end(), input_drf_option.begin(),  
                        [](unsigned char c){ return std::tolower(c); });

        if (input_drf_option == "equality") {
            return m_DetermineResourceAllocationByEquality(query_data, query_k);
        } else if (input_drf_option == "random") {
            return m_DetermineResourceAllocationByRandom(query_data, query_k);
        } else if (input_drf_option == "public") {
            return m_DetermineResourceAllocationByPublic(query_data, query_k);
        } else {
            std::cerr << "Unknown option for dominant resource allocation: " << drf_option 
                        << ", possible options are [equality | random | public]" << std::endl;
            std::exit(EXIT_FAILURE);            
        }
    }

    std::vector<size_t> m_DetermineResourceAllocationByPublic(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<size_t> allocation(m_silo_num, query_k);
        return allocation;
    }

    std::vector<size_t> m_DetermineResourceAllocationByGreedy(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option, const size_t& total_answer, const size_t& round_id, const std::set<std::pair<VectorDimensionType,int>>& minQ) {
        std::vector<size_t> allocation(m_silo_num, (size_t)0);
        if (m_silo_num==0 || query_k==0 || total_answer>=query_k) return allocation;

        if (round_id == 0) {
            std::fill(allocation.begin(), allocation.end(), (size_t)1);
            return allocation;
        }

        if (minQ.empty()) {
            std::cerr << "Tria's resource allocation: minQ shouldn't be empty" << std::endl;
            std::exit(EXIT_FAILURE);  
        }

        int silo_id = minQ.begin()->second;
        allocation[silo_id] = 1;

        return allocation;
    }

    std::vector<size_t> m_DetermineResourceAllocationByEquality(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<size_t> allocation(m_silo_num, (size_t)0);
        if (m_silo_num==0 || query_k==0) return allocation;

        // each data silo randomly gets equal amount
        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            allocation[silo_id] = query_k / m_silo_num;
        }

        // as for the remaining amounts, they should be randomly allocated to each data silo with equal probability
        const size_t remain_k = query_k % m_silo_num;
        if (remain_k > 0) {
            std::random_device rd;  
            std::mt19937 g(rd());  
            std::uniform_int_distribution<> dis(0, m_silo_num-1);  
            
            size_t sample_num = 0;
            while (sample_num++ < remain_k) {  
                size_t silo_id = dis(g);  
                ++allocation[silo_id];
            }
        }

        return allocation;
    }

    std::vector<size_t> m_DetermineResourceAllocationByRandom(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<size_t> allocation(m_silo_num, (size_t)0);
        if (m_silo_num==0 || query_k==0) return allocation;

        std::random_device rd;  
        std::mt19937 g(rd());  
        std::uniform_int_distribution<> dis(0, m_silo_num-1);  
        
        size_t sample_num = 0;
        while (sample_num++ < query_k) {  
            size_t silo_id = dis(g);  
            ++allocation[silo_id];
        }

        return allocation;
    }

    void m_CheckDrfOption(const std::string& drf_option) {
        std::string input_drf_option(drf_option);
        std::transform(input_drf_option.begin(), input_drf_option.end(), input_drf_option.begin(),  
                        [](unsigned char c){ return std::tolower(c); });

        if (input_drf_option == "equality") {
            std::cout << "Using ``equality`` dominant resource allocation" << std::endl;
        } else if (input_drf_option == "random") {
            std::cout << "Using ``random`` dominant resource allocation" << std::endl;
        } else if (input_drf_option == "public") {
            std::cout << "Using ``public`` as plaintext baseline" << std::endl;
        } else if (input_drf_option == "greedy") {
            std::cout << "Using ``greedy`` as our greedy algorithm" << std::endl;
        } else if (input_drf_option == "greedyx") {
            std::cout << "Using ``greedyX`` as our fast-greedy algorithm" << std::endl;
        } else {
            std::cerr << "Found unknown option for dominant resource allocation: " << drf_option 
                        << ", possible options are [equality | random | public | greedy | greedyX]" << std::endl;
            std::exit(EXIT_FAILURE);            
        }
    }

    void m_CreateSiloReceiver() {
        m_silo_num = m_silo_ipaddr_list.size();
        m_silo_receiver_list.resize(m_silo_num);
        grpc::ChannelArguments args;  
        args.SetInt(GRPC_ARG_MAX_SEND_MESSAGE_LENGTH, INT_MAX);  
        args.SetInt(GRPC_ARG_MAX_RECEIVE_MESSAGE_LENGTH, INT_MAX); 
        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            std::string silo_ipaddr = m_silo_ipaddr_list[silo_id];
            std::cout << "Server is connecting Silo #(" << std::to_string(silo_id+1) << ") on IP address " << silo_ipaddr << std::endl;

            std::shared_ptr<grpc::Channel> channel = grpc::CreateCustomChannel(silo_ipaddr, grpc::InsecureChannelCredentials(), args);
            m_silo_receiver_list[silo_id] = std::make_shared<SiloReceiver>(channel, silo_id+1, silo_ipaddr);
        }
    }

    std::vector<VectorDataType> m_PickGlobalKnn(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<VectorDataType> answer;
        if (std::is_same<Distance, InnerProductDistance>::value) {
            answer = m_PickGlobalKnn_ForMIPS(query_data, query_k);
        } else {
            answer = m_PickGlobalKnn_ByDefault(query_data, query_k);
        }
        return answer;
    }

    std::vector<VectorDataType> m_PickGlobalKnn_ByDefault(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<VectorDataType> query_answer;
        std::vector<VectorDataType> cand_answer;
        std::priority_queue<std::pair<VectorDimensionType, size_t>, std::vector<std::pair<VectorDimensionType, size_t>>, std::less<std::pair<VectorDimensionType, size_t>>> min_heap;  
        size_t min_heap_size = 0;
        Distance dist_function;

        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
            cand_answer.insert(cand_answer.end(), partial_knn_list.begin(), partial_knn_list.end());
        }
        for (size_t idx=0, cand_size=cand_answer.size(); idx<cand_size; ++idx) {  
            const auto& vector_data = cand_answer[idx];
            VectorDimensionType dist = dist_function(query_data, vector_data);
    
            if (min_heap_size < query_k) {
                ++min_heap_size;
                min_heap.push(std::make_pair(dist, idx));
            } else if (dist < min_heap.top().first) {  
                min_heap.pop();
                min_heap.push(std::make_pair(dist, idx));
            }  
        }

        while (!min_heap.empty()) {  
            size_t idx = min_heap.top().second; 
            query_answer.emplace_back(cand_answer[idx]);
            min_heap.pop();
        }  
    
        std::reverse(query_answer.begin(), query_answer.end());  

        return query_answer;
    }

    // specific for maximum inner product search
    std::vector<VectorDataType> m_PickGlobalKnn_ForMIPS(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<VectorDataType> query_answer;
        std::vector<VectorDataType> cand_answer;
        std::priority_queue<std::pair<VectorDimensionType, size_t>, std::vector<std::pair<VectorDimensionType, size_t>>, std::greater<std::pair<VectorDimensionType, size_t>>> min_heap;  
        size_t min_heap_size = 0;
        Distance dist_function;

        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[silo_id]->GetPartialKnn();
            cand_answer.insert(cand_answer.end(), partial_knn_list.begin(), partial_knn_list.end());
        }
        for (size_t idx=0, cand_size=cand_answer.size(); idx<cand_size; ++idx) {  
            const auto& vector_data = cand_answer[idx];
            VectorDimensionType dist = dist_function(query_data, vector_data);  
     
            if (min_heap_size < query_k) {
                ++min_heap_size;
                min_heap.push(std::make_pair(dist, idx));
            } else if (dist > min_heap.top().first) {  
                min_heap.pop();
                min_heap.push(std::make_pair(dist, idx));
            }  
        }
     
        while (!min_heap.empty()) {  
            size_t idx = min_heap.top().second; 
            query_answer.emplace_back(cand_answer[idx]);
            min_heap.pop();
        }  
    
        std::reverse(query_answer.begin(), query_answer.end());  

        return query_answer;
    }

    bool m_IsNotPlaintext(const std::string& drf_option) {
        std::string input_drf_option(drf_option);
        std::transform(input_drf_option.begin(), input_drf_option.end(), input_drf_option.begin(),  
                        [](unsigned char c){ return std::tolower(c); });
        return input_drf_option != "public";
    }

    bool m_IsBaselineMode(const std::string& drf_option) {
        std::string input_drf_option(drf_option);
        std::transform(input_drf_option.begin(), input_drf_option.end(), input_drf_option.begin(),  
                        [](unsigned char c){ return std::tolower(c); });
        return input_drf_option.find("greedy") == std::string::npos;
    }

    void m_InitBenchLogger() {
        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            m_silo_receiver_list[silo_id]->InitBenchLogger();
        }    
        m_logger.Init();    
    }

    std::vector<std::shared_ptr<SiloReceiver>> m_silo_receiver_list;
    std::string m_silo_ip_filename;
    std::string m_query_filename;
    std::vector<std::string> m_silo_ipaddr_list;
    std::vector<VectorDataType> m_query_data_list;
    std::vector<size_t> m_query_k_list;
    BenchLogger m_logger;
    int m_silo_num;
    const int m_thread_num;
};

std::unique_ptr<FedSqlServer> fed_sqlserver_ptr = nullptr;

void RunService(const std::string& silo_ip_filename, const std::string& query_filename, const std::string& output_filename, const std::string& truth_filename, const std::string& drf_option, const int& query_k) {
    fed_sqlserver_ptr = std::make_unique<FedSqlServer>(silo_ip_filename, query_filename, output_filename, truth_filename, drf_option, query_k);

    std::string log_info = fed_sqlserver_ptr->to_string();
    std::cout << log_info;
    std::cout.flush();
}    

// Ensure the log file is output, when the program is terminated.
void SignalHandler(int signal) {
    if (fed_sqlserver_ptr != nullptr) {
        std::string log_info = fed_sqlserver_ptr->to_string();
        std::cout << log_info;
        std::cout.flush();
    }
    quick_exit(0);
}

void ResetSignalHandler() {
    signal(SIGINT, SignalHandler);
    signal(SIGQUIT, SignalHandler);
    signal(SIGTERM, SignalHandler);
    signal(SIGKILL, SignalHandler);
}

int main(int argc, char** argv) {
    // Expect only arg: --silo-ip=../configuration/ip.txt --query-path=../query/query.txt --output-path=../query/answer.txt
    std::string query_filename, output_filename, truth_filename, silo_ip_filename;
    std::string drf_option("equality");
    int query_k;

    try { 
        bpo::options_description option_description("Required options");
        option_description.add_options()
            ("help", "produce help message")
            ("silo-ip", bpo::value<std::string>(), "Data silo's IP configuration file path")
            ("query-path", bpo::value<std::string>(), "Query file path")
            ("output-path", bpo::value<std::string>(), "Answer file path")
            ("truth-path", bpo::value<std::string>(), "Ground truth file path")
            ("option", bpo::value<std::string>(), "Option for k allocation: [equality | random | public | greedy | greedyX]")
            ("query-k", bpo::value<int>(&query_k)->default_value(128), "The K of KNN")
        ;

        bpo::variables_map variable_map;
        bpo::store(bpo::parse_command_line(argc, argv, option_description), variable_map);
        bpo::notify(variable_map);    

        if (variable_map.count("help")) {
            std::cout << option_description << std::endl;
            return 0;
        }

        bool options_all_set = true;

        if (variable_map.count("silo-ip")) {
            silo_ip_filename = variable_map["silo-ip"].as<std::string>();
            std::cout << "Data silo's IP configuration file path was set to " << silo_ip_filename << "\n";
        } else {
            std::cout << "Data silo's IP configuration file path was not set" << "\n";
            options_all_set = false;
        }

        if (variable_map.count("query-path")) {
            query_filename = variable_map["query-path"].as<std::string>();
            std::cout << "Server's received query file path was set to " << query_filename << "\n";
        } else {
            std::cout << "Server's received query file path was not set" << "\n";
            options_all_set = false;
        }

        if (variable_map.count("output-path")) {
            output_filename = variable_map["output-path"].as<std::string>();
            std::cout << "Server's answer output file path was set to " << output_filename << "\n";
        } else {
            std::cout << "Server's answer output file path was not set" << "\n";
        }

        if (variable_map.count("truth-path")) {
            truth_filename = variable_map["truth-path"].as<std::string>();
            std::cout << "Server's ground truth file path was set to " << truth_filename << "\n";
        } else {
            std::cout << "Server's ground truth file path was not set" << "\n";
        }

        if (variable_map.count("option")) {
            drf_option = variable_map["option"].as<std::string>();
            std::cout << "Server's k allocation strategy was set to " << drf_option << "\n";
        } else {
            std::cout << "Server's k allocation strategy was not set" << "\n";
        }

        if (variable_map.count("query-k")) {
            query_k = variable_map["query-k"].as<int>();
            std::cout << "The K of KNN was set to " << query_k << "\n";
        } else {
            std::cout << "The K of KNN was not set" << "\n";
        }

        if (false == options_all_set) {
            throw std::invalid_argument("Some options were not properly set");
            std::cout.flush();
            std::exit(EXIT_FAILURE);
        }

    } catch (std::exception& e) {  
        std::cerr << "Error: " << e.what() << "\n";  
        std::exit(EXIT_FAILURE);
    }

    ResetSignalHandler();

    RunService(silo_ip_filename, query_filename, output_filename, truth_filename, drf_option, query_k);

    return 0;
}

