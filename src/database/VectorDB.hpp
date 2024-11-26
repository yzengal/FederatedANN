/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.24
*/
#ifndef DATABASE_VECTORDB_HPP
#define DATABASE_VECTORDB_HPP

#include <iostream>  
#include <vector>
#include <memory>
#include <cstdlib>
#include <sstream>
#include <queue>
#include <utility>
#include <omp.h>

#include "../utils/DataType.hpp"
#include "../utils/MetricType.hpp"
#include "../utils/BenchLogger.hpp"
#include "../utils/File_IO.h"
#include "../index/BaseIndex.hpp"
#include "../index/LinearScanIndex.hpp"
#include "../index/HnswIndex.hpp"

template <typename Distance = EuclideanDistance>
class VectorDatabase {  
public:  
    explicit VectorDatabase(const std::string& file_name="") : m_construct_time(0), m_data_size(0) {
        m_logger.Init();

        if (!file_name.empty()) {
            ImportData(file_name);
        }

        omp_set_num_threads(1);
    }

    ~VectorDatabase() {

    }

    // Load the data from .csv or any other files
    void ImportData(const std::string& file_name) {
        m_logger.SetStartTimer();

        ReadVectorData(file_name, m_vector_data);
        m_data_size = m_vector_data.size();

        m_logger.SetEndTimer();
        double time_cost = m_logger.GetDurationTime() / 1000.0;
        std::cout << "Import data takes " << std::fixed << std::setprecision(3) << time_cost << " [s]" << std::endl;
    }

    // Build the index
    void ConstructIndex(const std::string& index_option) {
        m_logger.SetStartTimer();

        if (m_vector_index != nullptr)
            m_vector_index.reset();
        if (index_option.find("LinearScan") != std::string::npos) {
            m_vector_index = std::make_unique<LinearScanIndex<Distance>>(m_vector_data);
        } else if (index_option.find("HNSW") != std::string::npos) {
            m_vector_index = std::make_unique<HnswIndex<Distance>>(m_vector_data);
        } else {
            PrintLine(__LINE__);  
            std::cerr << "Unknown index option: " << index_option << ", valid options: [LinearScan | HNSW]" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        m_logger.SetEndTimer();
        double time_cost = m_logger.GetDurationTime() / 1000.0;
        m_construct_time = time_cost;
        std::cout << "Construct index takes " << std::fixed << std::setprecision(3) << time_cost << " [s]" << std::endl;
    } 

    // Query the input embedding and retrieve the k nearest neighbors
    std::vector<VectorDataType> KnnQuery(const VectorDataType& query_data, const size_t& query_k) {
        m_logger.SetStartTimer();

        std::vector<VectorDataType> query_answer;

        if (m_vector_index != nullptr) {
            query_answer = m_vector_index->KnnQuery(query_data, query_k);
        } else {
            if (std::is_same<Distance, InnerProductDistance>::value) {
                query_answer = m_KnnQueryByDefault_ForMIPS(query_data, query_k);
            } else {
                query_answer = m_KnnQueryByDefault(query_data, query_k);
            }
        }

        m_logger.SetEndTimer();
        m_logger.LogOneQuery();

        return query_answer;
    }

    // Query the input embedding and retrieve the k nearest neighbors in a batch
    std::vector<std::vector<VectorDataType>> BatchKnnQuery(const std::vector<VectorDataType>& query_data_list, const size_t& query_k) {
        m_logger.SetStartTimer();

        std::vector<std::vector<VectorDataType>> query_answer_list;
        if (m_vector_index != nullptr) {
            query_answer_list = m_vector_index->BatchKnnQuery(query_data_list, query_k);
        } else {
            if (std::is_same<Distance, InnerProductDistance>::value) {
                for (const auto& query_data : query_data_list) {
                    std::vector<VectorDataType> query_answer = m_KnnQueryByDefault_ForMIPS(query_data, query_k);
                    query_answer_list.emplace_back(query_answer);
                }
            } else {
                for (const auto& query_data : query_data_list) {
                    std::vector<VectorDataType> query_answer = m_KnnQueryByDefault(query_data, query_k);
                    query_answer_list.emplace_back(query_answer);
                }
            }
        }

        m_logger.SetEndTimer();
        m_logger.LogOneQuery();

        return query_answer_list;
    }

    // Output the current index size
    size_t IndexSize() const {
        if (m_vector_index == nullptr)
            return 0;
        return m_vector_index->IndexSize();
    }

    // Output the current data size
    size_t DataSize() const {
        return m_data_size;
    }

    double IndexTime() const {
        return m_construct_time;
    }

    std::string to_string() const {
        std::stringstream ss;

        ss << "data_size: " << DataSize() << "\n";
        ss << "construct_time: " << IndexTime() << " [s]\n";
        ss << "index_size: " << IndexSize() << " [kb]\n";
        ss << m_logger.to_string();

        return ss.str();
    }

private:
    std::vector<VectorDataType> m_KnnQueryByDefault(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<VectorDataType> query_answer;
        std::priority_queue<std::pair<VectorDimensionType, size_t>, std::vector<std::pair<VectorDimensionType, size_t>>, std::less<std::pair<VectorDimensionType, size_t>>> min_heap;  
        size_t min_heap_size = 0;
        Distance dist_function;

        for (size_t idx=0; idx<m_data_size; ++idx) {  
            const auto& vector_data = m_vector_data[idx];
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
            query_answer.emplace_back(m_vector_data[idx]);
            min_heap.pop();
        }  
    
        std::reverse(query_answer.begin(), query_answer.end());  

        return query_answer;
    }

    // specific for maximum inner product search
    std::vector<VectorDataType> m_KnnQueryByDefault_ForMIPS(const VectorDataType& query_data, const size_t& query_k) {
        std::vector<VectorDataType> query_answer;
        std::priority_queue<std::pair<VectorDimensionType, size_t>, std::vector<std::pair<VectorDimensionType, size_t>>, std::greater<std::pair<VectorDimensionType, size_t>>> min_heap;  
        size_t min_heap_size = 0;
        Distance dist_function;

        for (size_t idx=0; idx<m_data_size; ++idx) {  
            const auto& vector_data = m_vector_data[idx];
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
            query_answer.emplace_back(m_vector_data[idx]);
            min_heap.pop();
        }  
    
        std::reverse(query_answer.begin(), query_answer.end());  

        return query_answer;
    }

    std::unique_ptr<BaseIndex<Distance>> m_vector_index;
    std::vector<VectorDataType> m_vector_data;
    size_t m_data_size;
    double m_construct_time;
    BenchLogger m_logger;
};

#endif // DATABASE_VECTORDB_HPP