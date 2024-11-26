/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.29
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
#include <cstdlib>
#include <exception>
#include <signal.h>
#include <unistd.h>
#include <omp.h>

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

BenchLogger m_logger;

class FedVectorDBImpl final : public FedSqlService::Service {
public:
    explicit FedVectorDBImpl(const int silo_id, const std::string& silo_ipaddr, const std::string& data_filename, const std::string& index_type) {
        m_silo_ptr = std::make_unique<FileSiloConnector<Distance>>(silo_id, silo_ipaddr);
        
        m_silo_ptr->ImportData(data_filename);
        size_t data_size = m_silo_ptr->DataSize();
        std::cout << "Data silo holds " << std::to_string(data_size) << " vector data" << std::endl;

        m_logger.SetStartTimer();
        if(m_silo_ptr->IndexExists(data_filename, index_type)) {
            m_silo_ptr->LoadIndex(data_filename, index_type);
            std::cout << "Index loaded from disk" << std::endl;
        } else {
            m_silo_ptr->ConstructIndex(index_type);
        }
        m_logger.SetEndTimer();
        double construct_time = m_logger.GetDurationTime();
        std::cout << std::fixed << std::setprecision(6) << "Index construct time = " << construct_time/1000.0 << " [s]" << std::endl;
        size_t index_size = m_silo_ptr->IndexSize();
        std::cout << "Index size: " << index_size/1024 << " [KB]" << std::endl;
    }

    void SaveIndex(const std::string& file_path, const std::string& index_option) {
        m_silo_ptr->SaveIndex(file_path, index_option);
    }

    Status ExecuteKnnQuery(ServerContext* context,
                        const KnnSql* request,
                        Empty* empty_response) override {
        m_logger.SetStartTimer();

        #ifdef LOCAL_DEBUG
        std::cout << "Silo will perform ``ExecuteKnnQuery``" << std::endl;
        #endif

        const size_t dim = request->data_size();
        VectorDataType query_data(dim, request->vid());
        for (int i=0; i<dim; ++i) {
            query_data[i] = request->data(i);
        }
        size_t query_k = request->k();

        #ifdef LOCAL_DEBUG
        std::cout << "k = " << query_k << ", query_data = " << query_data.to_string() << std::endl;
        #endif

        m_local_knn_list = m_silo_ptr->KnnQuery(query_data, query_k);
        m_local_knn_head = 0;

        double grpc_comm = request->ByteSizeLong() + empty_response->ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);

        m_logger.SetEndTimer();
        m_logger.LogAddTime();

        #ifdef LOCAL_DEBUG
        std::cout << "Silo local knn:";
        for (auto vector_data : m_local_knn_list) {
            std::cout << " " << vector_data.vid;
        }
        std::cout << std::endl;
        #endif

        return Status::OK;
    }

    Status GetKnnCandidate(ServerContext* context,
                        const CandSql* request,
                        ServerWriter<VectorDataCand>* writer) override {
        m_logger.SetStartTimer();
        double grpc_comm = request->ByteSizeLong();

        const size_t cand_num = m_local_knn_list.size();
        size_t retrieve_num = request->k();

        if (m_local_knn_head+retrieve_num > cand_num) {
            throw std::invalid_argument("vector data silo has insufficient candidates for knn");
            std::exit(EXIT_FAILURE);
        }

        for (size_t i=0; i<retrieve_num; ++i) {
            VectorDataType& local_knn = m_local_knn_list[m_local_knn_head++];
            VectorDataCand vector_data_cand;
            vector_data_cand.set_vid(local_knn.vid);
            for (auto d : local_knn.data) {
                vector_data_cand.add_data(d);
            }
            writer->Write(vector_data_cand);
            grpc_comm += vector_data_cand.ByteSizeLong();
        }

        m_logger.LogAddComm(grpc_comm);
        m_logger.SetEndTimer();
        m_logger.LogAddTime();

        return Status::OK;
    }

    Status ClearKnnCandidate(ServerContext* context,
                        const Empty* request,
                        Empty* empty_response) override {

        m_logger.SetStartTimer();

        m_local_knn_list.clear();
        m_local_knn_head = 0;

        double grpc_comm = request->ByteSizeLong() + empty_response->ByteSizeLong();
        m_logger.LogAddComm(grpc_comm);

        m_logger.SetEndTimer();
        m_logger.LogOneQuery();

        return Status::OK;
    }

    std::string to_string() const {
        std::stringstream ss;

        ss << "-------------- Data Silo #(" << m_silo_ptr->GetSiloId() << ") Log --------------\n";
        ss << m_logger.to_string();

        return ss.str();
    }

private:
    std::unique_ptr<BaseSiloConnector<Distance>> m_silo_ptr;
    std::vector<VectorDataType> m_local_knn_list; // the most similar one appears at the head, the least similar one appears at the tail
    size_t m_local_knn_head;                      // the head pointer of the vector m_local_knn_list
    BenchLogger m_logger;
};

std::unique_ptr<FedVectorDBImpl> fed_vectordb_ptr = nullptr;

void RunSilo(const int silo_id, const std::string& silo_ipaddr, const std::string& table_filename, const std::string& index_type) {
    fed_vectordb_ptr = std::make_unique<FedVectorDBImpl>(silo_id, silo_ipaddr, table_filename, index_type);

    ServerBuilder builder;
    builder.AddListeningPort(silo_ipaddr, grpc::InsecureServerCredentials());
    builder.RegisterService(fed_vectordb_ptr.get());
    builder.SetMaxSendMessageSize(INT_MAX);
    builder.SetMaxReceiveMessageSize(INT_MAX);
    std::unique_ptr<Server> server(builder.BuildAndStart());
    fed_vectordb_ptr->SaveIndex(table_filename, index_type);
    std::cout << "Data Silo #(" << silo_id << ") listening on " << silo_ipaddr << std::endl;
    
    server->Wait();

    std::string log_info = fed_vectordb_ptr->to_string();
    std::cout << log_info;
    std::cout.flush();
}

// Ensure the log file is output, when the program is terminated.
void SignalHandler(int signal) {
    if (fed_vectordb_ptr != nullptr) {
        std::string log_info = fed_vectordb_ptr->to_string();
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
    omp_set_num_threads(4);
    
    // Expect the following args: --ip=0.0.0.0 --port=50051 --data-path=../../data/data_01.txt --silo_id=1
    int silo_port, silo_id;
    std::string silo_ip, silo_ipaddr, data_filename, index_type;
    
    try { 
        bpo::options_description option_description("Required options");
        option_description.add_options()
            ("help", "produce help message")
            ("id", bpo::value<int>(&silo_id)->default_value(0), "Data silo's ID")
            ("ip", bpo::value<std::string>(), "Data silo's IP address")
            ("port", bpo::value<int>(&silo_port), "Data silo's IP port")
            ("data-path", bpo::value<std::string>(), "Data file path")
            ("index-type", bpo::value<std::string>(), "Type of index type: [Linear | HNSW | LSH | PQ | IVF | IVFPQ]")
        ;

        bpo::variables_map variable_map;
        bpo::store(bpo::parse_command_line(argc, argv, option_description), variable_map);
        bpo::notify(variable_map);    

        if (variable_map.count("help")) {
            std::cout << option_description << std::endl;
            return 0;
        }

        bool options_all_set = true;

        std::cout << "Data silo's ID is " << silo_id << "\n";

        if (variable_map.count("ip")) {
            silo_ip = variable_map["ip"].as<std::string>();
            std::cout << "Data silo's IP address was set to " << silo_ip << "\n";
        } else {
            std::cout << "Data silo's IP address was not set" << "\n";
            options_all_set = false;
        }

        if (variable_map.count("port")) {
            silo_port = variable_map["port"].as<int>();
            std::cout << "Data silo's IP port was set to " << silo_port << "\n";
        } else {
            std::cout << "Data silo's IP port was not set" << "\n";
            options_all_set = false;
        }

        if (variable_map.count("data-path")) {
            data_filename = variable_map["data-path"].as<std::string>();
            std::cout << "Data silo's data file path was set to " << data_filename << "\n";
        } else {
            std::cout << "Data silo's data file path was not set" << "\n";
            options_all_set = false;
        }

        if (variable_map.count("index-type")) {
            index_type = variable_map["index-type"].as<std::string>();
            std::cout << "Data silo's index type was set to " << data_filename << "\n";
        } else {
            std::cout << "Data silo's index type was not set" << "\n";
            options_all_set = false;
        }

        if (false == options_all_set) {
            throw std::invalid_argument("Some options were not properly set");
            std::cout.flush();
            std::exit(EXIT_FAILURE);
        }

        silo_ipaddr = silo_ip + std::string(":") + std::to_string(silo_port);
    } catch (std::exception& e) {  
        std::cerr << "Error: " << e.what() << "\n";  
        std::exit(EXIT_FAILURE);
    }

    ResetSignalHandler();

    RunSilo(silo_id, silo_ipaddr, data_filename, index_type);

    return 0;
}

