/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.23
*/
#ifndef INDEX_BASE_INDEX_HPP
#define INDEX_BASE_INDEX_HPP

#include <iostream>
#include <vector>

#include "../utils/DataType.hpp"
#include "../utils/MetricType.hpp"

template <typename Distance = EuclideanDistance>
class BaseIndex {  
public:
    virtual ~BaseIndex() = default;

    // Construct the index
    virtual void Construct(const std::vector<VectorDataType>& data_list) = 0;

    // Load the index
    virtual void Load(const std::string& index_path, const std::vector<VectorDataType>& data_list) = 0;

    // Query the input embedding and retrieve the k nearest neighbors
    virtual std::vector<VectorDataType> KnnQuery(const VectorDataType& query_data, const size_t& query_k) = 0;

    // Query the input embedding and retrieve the k nearest neighbors along with their distances
    virtual void KnnQuery(const VectorDataType& query_data, const size_t& query_k, std::vector<VectorDataType>& answer_object, std::vector<VectorDimensionType>& answer_dist) = 0;

    // Output the current index size
    virtual size_t IndexSize() const = 0;

    // Batch knn query
    virtual std::vector<std::vector<VectorDataType>> BatchKnnQuery(const std::vector<VectorDataType>& query_data_list, const size_t& query_k) {
        std::vector<std::vector<VectorDataType>> batch_answer_list;

        for (const auto& query_data : query_data_list) {
            std::vector<VectorDataType> answer_list = KnnQuery(query_data, query_k);
            batch_answer_list.emplace_back(answer_list);
        }

        return batch_answer_list;
    }

    // Output the current data size
    virtual size_t DataSize() const = 0;

    virtual void SaveIndex(const std::string& index_path) const = 0;
};

#endif // INDEX_BASE_INDEX_HPP