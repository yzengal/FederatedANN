/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.11.05
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

#include "FedVectorDBO.grpc.pb.h"

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

using FedVectorDBO::KnnSql;
using FedVectorDBO::CandSql;
using FedVectorDBO::MpcComm;
using FedVectorDBO::MinMpcBobList;
using FedVectorDBO::MinMpcAlicePair;
using FedVectorDBO::MinMpcSecretShares;
using FedVectorDBO::VectorDataCand;
using FedVectorDBO::FedSqlService;

using Distance = EuclideanSquareDistance;


class SiloReceiver {
public:
    SiloReceiver(std::shared_ptr<grpc::Channel> channel, const int silo_id, const std::string& silo_ipaddr) 
        : m_silo_stub(FedSqlService::NewStub(channel)), m_silo_id(silo_id), m_silo_ipaddr(silo_ipaddr) {
        
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

        Status status = m_silo_stub->ExecuteKnnQuery(&context, query, &response); 
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
            m_silo_stub->GetKnnCandidate(&context, query));

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

        Status status = m_silo_stub->ClearKnnCandidate(&context, request, &response); 
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

    void StartMinMpc(int bob_silo_id) {
        ClientContext context;
        MinMpcBobList request;
        Empty response;

        request.set_bob(bob_silo_id); // negative means all silos should be regarded as Bob
        Status status = m_silo_stub->StartMinMpc(&context, request, &response); 
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Start the protocol for Min-MPC from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }

        float grpc_comm = request.ByteSizeLong() + response.ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);
    }

    void CollectMinMpc(int bob_silo_id) {
        ClientContext context;
        MinMpcBobList request;
        MinMpcSecretShares response;

        request.set_bob(bob_silo_id); // negative means all silos should be regarded as Bob
        Status status = m_silo_stub->CollectMinMpc(&context, request, &response); 
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Get secret share list from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }

        if (bob_silo_id < 0) {
            m_secret_share_list.clear();
            int secret_share_num = response.r_size();
            for (int i=0; i<secret_share_num; ++i) {
                m_secret_share_list.emplace_back(response.r(i));
            }
        } else {
            assert(m_silo_id != bob_silo_id);
            int secret_share_num = response.r_size();
            assert(secret_share_num == 2);
            m_secret_share_list[m_silo_id] = response.r(0);
            m_secret_share_list[bob_silo_id] = response.r(1);
        }

        #ifdef LOCAL_DEBUG
        std::cout << "Server receives secret share from silo " << m_silo_id << ": [";
        for (int i=0; i<m_secret_share_list.size(); ++i) {
            if (i > 0)
                std::cout << ", ";
            std::cout << m_secret_share_list[i];
        }
        std::cout << "]" << std::endl;
        #endif

        float grpc_comm = request.ByteSizeLong() + response.ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);
    }

    std::vector<int> GetSecretShareList() const {
        return m_secret_share_list;
    }

    int GetSecretShare(int idx) const {
        return m_secret_share_list[idx];
    }

    std::vector<VectorDataType> GetPartialKnn() const {
        return m_partial_knn;
    }

    double GetMpcComm() const {
        return m_mpc_comm;
    }

    void CollectMpcComm() {
        double mpc_comm = 0.0;

        ClientContext context;
        Empty request;
        MpcComm response;

        Status status = m_silo_stub->GetMpcComm(&context, request, &response); 
        if (!status.ok()) {
            std::cerr << "RPC failed: " << status.error_message() << std::endl;
            std::string error_message;
            error_message = std::string("Get MPC communication from data silo #(") + std::to_string(m_silo_id) + std::string(") failed");
            throw std::invalid_argument(error_message);  
            std::exit(EXIT_FAILURE);
        }
        m_mpc_comm = response.comm();
    }

    double GetQueryComm() const {
        return m_logger.GetQueryComm();
    }

    void InitBenchLogger() {
        m_logger.Init();
    }

    static void StartMinMpcThread(SiloReceiver* silo_receiver, int bob_silo_id) {  
        silo_receiver->StartMinMpc(bob_silo_id);  
    } 

    static void CollectMinMpcThread(SiloReceiver* silo_receiver, int bob_silo_id) {  
        silo_receiver->CollectMinMpc(bob_silo_id);  
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
    std::unique_ptr<FedSqlService::Stub> m_silo_stub;
    std::vector<VectorDataType> m_partial_knn;
    std::vector<int> m_secret_share_list;
    std::string m_silo_ipaddr;
    double m_mpc_comm;
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
        m_InitLocalCache();
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
        #ifdef LOCAL_DEBUG
        query_num = 1;
        #endif
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
        }

        // Step 3: Print the log information
        m_logger.Print();
    }

    std::vector<VectorDataType> ProcessOneFedKnnQuery(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option="equality") {
        std::vector<VectorDataType> answer;
        
        answer = m_ProcessFedKnnQuery_ByGreedyO(query_data, query_k, drf_option);
        
        return answer;
    }

    std::string to_string() const {
        std::stringstream ss;

        ss << "-------------- Service Log --------------\n";
        ss << m_logger.to_string();

        return ss.str();
    }

private:
    std::vector<VectorDataType> m_ProcessFedKnnQuery_ByGreedyO(const VectorDataType& query_data, const size_t& query_k, const std::string& drf_option) {
        m_logger.SetStartTimer();
        Distance dist_function;

        // Step 0: Intialize some variables
        std::vector<VectorDataType> answer;

        // Step 1: Perform local knn query at each silo
        m_PerformLocalKnn(query_data, query_k);

        // Step 2: Determine the resource allocation for each silo
        std::vector<VectorDimensionType> dist_list(m_silo_num);
        const int retrieve_step_size = 1;
        int prev_optimal_silo_id = -1 - query_data.vid;
        for (size_t rid=0, total_answer=0; total_answer<query_k; ++rid) {
            int optimal_silo_id = m_GetMostRelevantSilo(rid, prev_optimal_silo_id);

            SiloReceiver::RequestPartialKnnThread(m_silo_receiver_list[optimal_silo_id].get(), retrieve_step_size);
            std::vector<VectorDataType> partial_knn_list = m_silo_receiver_list[optimal_silo_id]->GetPartialKnn();
            answer.insert(answer.end(), partial_knn_list.begin(), partial_knn_list.end());
            total_answer += retrieve_step_size;

            prev_optimal_silo_id = optimal_silo_id;
        }

        // Step 4: Clear Knn candidates at each silo
        m_ClearKnnCand();

        // Step 5: Print the log information
        m_logger.SetEndTimer();
        double mpc_comm = 0.0;
        double query_comm = 0.0;
        for (int i=0; i<m_silo_num; ++i) {
            m_silo_receiver_list[i]->CollectMpcComm();
            double cur_mpc_comm = m_silo_receiver_list[i]->GetMpcComm();
            mpc_comm += cur_mpc_comm;
        }
        query_comm += mpc_comm;
        for (int i=0; i<m_silo_num; ++i) {
            query_comm += m_silo_receiver_list[i]->GetQueryComm();
        }
        query_comm -= m_logger.GetQueryComm();
        double query_time = m_logger.GetDurationTime();
        m_logger.LogOneQuery(query_comm);

        std::cout << std::fixed << std::setprecision(6) 
                    << "Query #(" << query_data.vid << "): runtime = " << query_time << " [ms], communication = " << query_comm/1024.0 << " [KB]" << std::endl;
    
        return answer;
    }

    int m_GetMostRelevantSilo(int rid, int prev_optimal_silo_id) {
        int optimal_silo_id = -1;

        // 1. enable all data silo as alice
        m_StartMinMpc(prev_optimal_silo_id);

        // 2. collect the secret shares from all silos
        m_CollectMinMpc(prev_optimal_silo_id);

        // 3. pick the minimum (if there're ties, randomly pick on)
        m_AssembleSecretShare(prev_optimal_silo_id);

        int mn = INT_MAX;
        for (int i=0; i<m_silo_num; ++i) {
            if (m_sum_secret_list[i] < mn) {
                mn = m_sum_secret_list[i];
                optimal_silo_id = i;
            }
        }

        #ifdef LOCAL_DEBUG
        std::cout << "Sum of secret shares: [";
        for (int i=0; i<m_silo_num; ++i) {
            if (i > 0)
                std::cout << ", ";
            std::cout << m_sum_secret_list[i];
        }
        std::cout << "]" << std::endl;
        std::cout << "Most relevant silo: " << optimal_silo_id << ", value = " << std::fixed << std::setprecision(6) << mn << std::endl;
        std::cout << "------------------------------------------------" << std::endl;
        #endif

        return optimal_silo_id;
    }

    void m_AssembleSecretShare(int prev_optimal_silo_id) {
        if (prev_optimal_silo_id < 0) {

            std::fill(m_sum_secret_list.begin(), m_sum_secret_list.end(), 0);
            for (int i=0; i<m_silo_num; ++i) {
                m_secret_share_list[i] = m_silo_receiver_list[i]->GetSecretShareList();
                for (int j=0; j<m_silo_num; ++j) {
                    m_sum_secret_list[j] += m_secret_share_list[i][j];
                }
            }

        } else {

            for (int i=0; i<m_silo_num; ++i) {
                if (i == prev_optimal_silo_id) {
                    std::vector<int> tmp_secret_share_list = m_silo_receiver_list[i]->GetSecretShareList();
                    for (int j=0; j<m_silo_num; ++j) {
                        m_sum_secret_list[j] += (tmp_secret_share_list[j] - m_secret_share_list[i][j]);
                        m_secret_share_list[i][j] = tmp_secret_share_list[j];
                    }
                                   
                } else {

                    // update the secret share when $i$th silo acts as alice
                    int alice_secret_share = m_silo_receiver_list[i]->GetSecretShare(i);
                    m_sum_secret_list[i] += (alice_secret_share - m_secret_share_list[i][i]);
                    m_secret_share_list[i][i] = alice_secret_share;

                    // update the secret share when $i$th silo acts as bob
                    int bob_secret_share = m_silo_receiver_list[i]->GetSecretShare(prev_optimal_silo_id);
                    m_sum_secret_list[prev_optimal_silo_id] += (bob_secret_share - m_secret_share_list[i][prev_optimal_silo_id]);
                    m_secret_share_list[i][prev_optimal_silo_id] = bob_secret_share;
                }
            }

        }

        #ifdef LOCAL_DEBUG
        std::cout << "Server maintains sum of secret shares: [";
        for (int i=0; i<m_sum_secret_list.size(); ++i) {
            if (i > 0)
                std::cout << ", ";
            std::cout << m_sum_secret_list[i];
        }
        std::cout << "]" << std::endl;
        #endif
    }

    void m_StartMinMpc(int prev_optimal_silo_id) {
        if (prev_optimal_silo_id < 0) {
            boost::asio::thread_pool pool(m_thread_num);
            const int bob_silo_id = prev_optimal_silo_id;

            for (int alice_silo_id=0; alice_silo_id<m_silo_num; ++alice_silo_id) {
                boost::asio::post(pool, std::bind(SiloReceiver::StartMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 
            }

            pool.join(); 

        } else {
            int alice_silo_id, bob_silo_id;

            boost::asio::thread_pool pool(m_thread_num);

            // 1. currnt silo (i.e., prev_optimal_silo_id) acts as Alice and changes all
            alice_silo_id = prev_optimal_silo_id;
            bob_silo_id = -2;
            boost::asio::post(pool, std::bind(SiloReceiver::StartMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 

            // 2. currnt silo (i.e., prev_optimal_silo_id) acts as Bob and is changed by Alice
            bob_silo_id = prev_optimal_silo_id;
            for (alice_silo_id=0; alice_silo_id<m_silo_num; ++alice_silo_id) {
                if (alice_silo_id == bob_silo_id)
                    continue;
                boost::asio::post(pool, std::bind(SiloReceiver::StartMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 
            }            
     
            pool.join(); 
        }
    }

    void m_CollectMinMpc(const int prev_optimal_silo_id) {
        if (prev_optimal_silo_id < 0) {
            boost::asio::thread_pool pool(m_thread_num);
            const int bob_silo_id = prev_optimal_silo_id;

            for (int alice_silo_id=0; alice_silo_id<m_silo_num; ++alice_silo_id) {
                boost::asio::post(pool, std::bind(SiloReceiver::CollectMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 
            }

            pool.join(); 
                   
        } else {
            int alice_silo_id, bob_silo_id;

            boost::asio::thread_pool pool(m_thread_num);

            // 1. currnt silo (i.e., prev_optimal_silo_id) acts as Alice and changes all
            alice_silo_id = prev_optimal_silo_id;
            bob_silo_id = -1;
            boost::asio::post(pool, std::bind(SiloReceiver::CollectMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 

            // 2. currnt silo (i.e., prev_optimal_silo_id) acts as Bob and is changed by Alice
            bob_silo_id = prev_optimal_silo_id;
            for (alice_silo_id=0; alice_silo_id<m_silo_num; ++alice_silo_id) {
                if (alice_silo_id == bob_silo_id)
                    continue;
                boost::asio::post(pool, std::bind(SiloReceiver::CollectMinMpcThread, m_silo_receiver_list[alice_silo_id].get(), bob_silo_id)); 
            }            
     
            pool.join(); 
        }       
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

    void m_CheckDrfOption(const std::string& drf_option) {
        std::string input_drf_option(drf_option);
 
        if (input_drf_option == "greedyO") {
            std::cout << "Using ``optimized greedy`` as our greedy algorithm for k is smaller than m" << std::endl;
        } else {
            std::cerr << "Found unknown option for dominant resource allocation: " << drf_option 
                        << ", possible options are [greedyO]" << std::endl;
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
            std::cout << "Server is connecting Silo #(" << std::to_string(silo_id) << ") on IP address " << silo_ipaddr << std::endl;

            std::shared_ptr<grpc::Channel> channel = grpc::CreateCustomChannel(silo_ipaddr, grpc::InsecureChannelCredentials(), args);
            m_silo_receiver_list[silo_id] = std::make_shared<SiloReceiver>(channel, silo_id, silo_ipaddr);
        }
    }

    void m_InitLocalCache() {
        m_sum_secret_list.resize(m_silo_num, 0);
        m_secret_share_list.resize(m_silo_num);
        std::fill(m_sum_secret_list.begin(), m_sum_secret_list.end(), 0);
        for (int i=0; i<m_silo_num; ++i) {
            m_secret_share_list[i].resize(m_silo_num);
            std::fill(m_secret_share_list[i].begin(), m_secret_share_list[i].end(), 0);
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

    void m_InitBenchLogger() {
        for (int silo_id=0; silo_id<m_silo_num; ++silo_id) {
            m_silo_receiver_list[silo_id]->InitBenchLogger();
        }    
        m_logger.Init();    
    }

    std::vector<std::vector<int>> m_secret_share_list;
    std::vector<int> m_sum_secret_list;
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
    std::string drf_option("greedyO");
    int query_k;

    try { 
        bpo::options_description option_description("Required options");
        option_description.add_options()
            ("help", "produce help message")
            ("silo-ip", bpo::value<std::string>(), "Data silo's IP configuration file path")
            ("query-path", bpo::value<std::string>(), "Query file path")
            ("output-path", bpo::value<std::string>(), "Answer file path")
            ("truth-path", bpo::value<std::string>(), "Ground truth file path")
            ("option", bpo::value<std::string>(), "Option for k allocation: [greedyO]")
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

