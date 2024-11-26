/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.23
*/
#ifndef INDEX_HNSW_INDEX_HPP
#define INDEX_HNSW_INDEX_HPP

#include <iostream>  
#include <vector>
#include <stdexcept>
#include <faiss/impl/HNSW.h>
#include <faiss/IndexHNSW.h>
#include <faiss/index_io.h>
#include <faiss/index_io.h>

#include "BaseIndex.hpp"
#include "../utils/MetricType.hpp"
#include "../utils/BenchLogger.hpp"

template <typename Distance = EuclideanDistance>
class HnswIndex final : public BaseIndex<Distance> { 
public:  
    HnswIndex() : m_hnsw(NULL), m_data(NULL), m_dim(0), m_data_size(0)  {

    } 

    HnswIndex(const std::vector<VectorDataType>& data_list) : HnswIndex() {
        Construct(data_list);
    }

    HnswIndex(const std::string& index_path, const std::vector<VectorDataType>& data_list) : HnswIndex() {
        Load(index_path, data_list);
    }
    
    ~HnswIndex() override {
        delete m_hnsw;
    }  
  
    void Construct(const std::vector<VectorDataType>& data_list) override {
        m_ClearHnsw();
        if (data_list.empty()) return ;

        m_data = data_list.data();
        m_dim = data_list[0].Dimension();
        m_data_size = data_list.size();

        if (std::is_same<Distance, EuclideanDistance>::value) {
            m_hnsw = new faiss::IndexHNSWFlat(m_dim, m_hnsw_M, faiss::MetricType::METRIC_L2);
        } else if (std::is_same<Distance, EuclideanSquareDistance>::value) {
            m_hnsw = new faiss::IndexHNSWFlat(m_dim, m_hnsw_M, faiss::MetricType::METRIC_L2);
        } else if (std::is_same<Distance, InnerProductDistance>::value) {
            m_hnsw = new faiss::IndexHNSWFlat(m_dim, m_hnsw_M, faiss::MetricType::METRIC_INNER_PRODUCT);
        } else {
            throw std::invalid_argument("Metric distance does not support, and invalid ones are [L2 | inner_product]");
            std::exit(EXIT_FAILURE);
        }
        m_hnsw->hnsw.efConstruction = 500;
        m_hnsw->hnsw.efSearch = 100;

        float* xb = new float[m_dim * m_data_size];

        for (size_t i=0; i<m_data_size; ++i) {
            for (size_t j=0; j<m_dim; ++j) {
                xb[i*m_dim+j] = data_list[i].data[j];
            }
        }

        m_hnsw->add(m_data_size, xb); // add vectors to the index

        delete[] xb;
    }

    void Load(const std::string& index_path, const std::vector<VectorDataType>& data_list) override {
        faiss::Index* loadedIndex = faiss::read_index(index_path.c_str());
        m_hnsw = dynamic_cast<faiss::IndexHNSWFlat*>(loadedIndex);
        m_data = data_list.data();
        m_dim = m_hnsw->d;
        m_data_size = m_hnsw->ntotal;
        std::cout << "Load hnsw index, dim = " << m_dim << ", data_size = " << m_data_size << std::endl;
    }
    
    size_t IndexSize() const override {
        size_t index_size = 0;

        // faiss::IndexHNSWFlat *m_hnsw;
        // VectorDataType* m_data;
        index_size += sizeof(void *) * 2;

        // int d;        ///< vector dimension
        // idx_t ntotal; ///< total nb of indexed vectors
        // bool verbose; ///< verbosity level
        // bool is_trained;
        // MetricType metric_type;
        // float metric_arg;
        index_size += sizeof(int);
        index_size += sizeof(faiss::idx_t);
        index_size += sizeof(bool)*2;
        index_size += sizeof(faiss::MetricType);
        index_size += sizeof(float);

        // the link strcuture
        // HNSW hnsw;
        index_size += m_hnsw->hnsw.assign_probas.size() * sizeof(double) + sizeof(size_t);
        index_size += m_hnsw->hnsw.cum_nneighbor_per_level.size() * sizeof(int) + sizeof(size_t);
        index_size += m_hnsw->hnsw.levels.size() * sizeof(int) + sizeof(size_t);
        index_size += m_hnsw->hnsw.offsets.size() * sizeof(size_t) + sizeof(size_t);
        index_size += m_hnsw->hnsw.neighbors.size() * sizeof(faiss::HNSW::storage_idx_t) + sizeof(size_t);
        index_size += sizeof(faiss::HNSW::storage_idx_t);
        index_size += sizeof(int)*3;
        index_size += sizeof(bool);
        index_size += sizeof(int);
        index_size += sizeof(bool);

        return index_size;
    }

    size_t DataSize() const override {
        return m_data_size;
    }

    std::vector<VectorDataType> KnnQuery(const VectorDataType& query_data, const size_t& query_k) override {
        std::vector<VectorDataType> answer_object;
        std::vector<VectorDimensionType> answer_dist;

        m_KnnQuery(query_data, query_k, answer_object, answer_dist);
        
        return answer_object;
    }

    void KnnQuery(const VectorDataType& query_data, const size_t& query_k, std::vector<VectorDataType>& answer_object, std::vector<VectorDimensionType>& answer_dist) override {
        m_KnnQuery(query_data, query_k, answer_object, answer_dist);
    }

    std::vector<std::vector<VectorDataType>> BatchKnnQuery(const std::vector<VectorDataType>& query_data_list, const size_t& query_k) override {
        std::vector<std::vector<VectorDataType>> answer_object_list;
        std::vector<std::vector<VectorDimensionType>> answer_dist_list;

        m_BatchKnnQuery(query_data_list, query_k, answer_object_list, answer_dist_list);

        return answer_object_list;
    }

    void SaveIndex(const std::string& index_path) const override {
        std::cout << "Saving HNSW index..." << std::endl;
        faiss::write_index(m_hnsw, index_path.c_str());
    }

private:
    void m_KnnQuery(const VectorDataType& query_data, const size_t& query_k, std::vector<VectorDataType>& answer_object, std::vector<VectorDimensionType>& answer_dist) {

        answer_object.clear();
        answer_dist.clear();

        faiss::idx_t* idx_list = new faiss::idx_t[query_k];
        float* dist_list = new float[query_k];
        float* xq = new float[m_dim];

        for (size_t j=0; j<m_dim; ++j) {
            xq[j] = query_data.data[j];
        }

        std::fill(idx_list, idx_list+query_k, -1);
        std::fill(dist_list, dist_list+query_k, 0.0);

        #ifdef LOCAL_DEBUG
        PrintLine(__LINE__);
        std::cout << "k = " << query_k << ", query_data = " << query_data.to_string() << std::endl;
        std::vector<float> tmp_float_list(xq, xq+m_dim);
        PrintLine(__LINE__);
        std::cout << "m_dim = " << m_dim << ", xq = ";
        PrintVector<float>(tmp_float_list);
        #endif

        m_hnsw->search(1, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i=0; i<query_k; ++i) {
                if (idx_list[i]>=0 && dist_list[i]>=0)
                    dist_list[i] = std::sqrt(dist_list[i]);
            }
        }

        for (size_t i=0; i<query_k; ++i) {
            size_t idx = idx_list[i];
            answer_object.emplace_back(m_data[idx]);
            answer_dist.emplace_back(dist_list[i]);
        }

        delete[] idx_list;
        delete[] dist_list;
        delete[] xq;
    }

    void m_BatchKnnQuery(const std::vector<VectorDataType>& query_data_list, const size_t& query_k, std::vector<std::vector<VectorDataType>>& answer_object_list, std::vector<std::vector<VectorDimensionType>>& answer_dist_list) {
        answer_object_list.clear();
        answer_dist_list.clear();

        const size_t nq = query_data_list.size();
        const size_t answer_size = query_k * nq;

        faiss::idx_t* idx_list = new faiss::idx_t[answer_size];
        float* dist_list = new float[answer_size];
        float* xq = new float[m_dim * nq];

        for (size_t i=0; i<nq; ++i) {
            for (size_t j=0; j<m_dim; ++j) {
                xq[i*m_dim+j] = query_data_list[i].data[j];
            }
        }

        std::fill(idx_list, idx_list+answer_size, -1);
        std::fill(dist_list, dist_list+answer_size, 0.0);

        m_hnsw->search(nq, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i=0; i<answer_size; ++i) {
                if (idx_list[i]>=0 && dist_list[i]>=0)
                    dist_list[i] = std::sqrt(dist_list[i]);
            }
        }

        for (size_t i=0; i<nq; ++i) {
            std::vector<VectorDataType> knn_list;
            std::vector<VectorDimensionType> knnd_list;
            for (size_t j=0; j<query_k; ++j) {
                size_t idx = idx_list[i*query_k+j];
                knn_list.emplace_back(m_data[idx]);
                knnd_list.emplace_back(dist_list[i*query_k+j]);
            }
            answer_object_list.emplace_back(knn_list);
            answer_dist_list.emplace_back(knnd_list);
        }

        delete[] idx_list;
        delete[] dist_list;
        delete[] xq;
    }

    void m_ClearHnsw() {
        if (m_hnsw != NULL)
            delete m_hnsw;
        m_hnsw = NULL;
        m_data = NULL;
    }

    faiss::IndexHNSWFlat *m_hnsw;
    const VectorDataType* m_data;
    int m_dim;
    size_t m_data_size;
    const int m_hnsw_M = 32;
};

#endif // INDEX_HNSW_INDEX_HPP