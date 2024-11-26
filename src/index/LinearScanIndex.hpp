/**
	@author:	Yuxiang Zeng
	@email: 	yxzeng@buaa.edu.cn
	@date:		2024.09.23
*/
#ifndef INDEX_LINEARSCAN_INDEX_HPP
#define INDEX_LINEARSCAN_INDEX_HPP

#include <iostream>
#include <vector>
#include <faiss/IndexFlat.h>
#include <faiss/index_io.h>
#include <faiss/index_io.h>

#include "BaseIndex.hpp"
#include "../utils/MetricType.hpp"

template<typename Distance = EuclideanDistance>
class LinearScanIndex final : public BaseIndex<Distance> {
public:
    LinearScanIndex() : m_data(NULL), m_dim(0), m_data_size(0) {

    }

    LinearScanIndex(const std::vector <VectorDataType> &data_list) : LinearScanIndex() {
        Construct(data_list);
    }

    LinearScanIndex(const std::string& index_path, const std::vector<VectorDataType>& data_list) : LinearScanIndex() {
        Load(index_path, data_list);
    }

    ~LinearScanIndex() override {
        delete m_index_flat;
    }

    void Construct(const std::vector<VectorDataType>& data_list) override {
        //m_ClearIndexFlat();
        if (data_list.empty()) return;

        m_data = data_list.data();
        m_dim = data_list[0].Dimension();
        m_data_size = data_list.size();

        if (std::is_same<Distance, EuclideanDistance>::value) {
            m_index_flat = new faiss::IndexFlatL2(m_dim);
        } else if (std::is_same<Distance, EuclideanSquareDistance>::value) {
            m_index_flat = new faiss::IndexFlatL2(m_dim);
        } else if (std::is_same<Distance, InnerProductDistance>::value) {
            m_index_flat = new faiss::IndexFlatIP(m_dim);
        } else {
            throw std::invalid_argument("Metric distance does not support, and invalid ones are [L2 | inner_product]");
            std::exit(EXIT_FAILURE);
        }

        float *xb = new float[m_dim * m_data_size];

        for (size_t i = 0; i < m_data_size; ++i) {
            for (size_t j = 0; j < m_dim; ++j) {
                xb[i * m_dim + j] = data_list[i].data[j];
            }
        }

        m_index_flat->add(m_data_size, xb); // add vectors to the index

        delete[] xb;
    }

    void Load(const std::string& index_path, const std::vector<VectorDataType>& data_list) override {
        faiss::Index* loadedIndex = faiss::read_index(index_path.c_str());
        m_index_flat = dynamic_cast<faiss::IndexFlat*>(loadedIndex);
        m_data = data_list.data();
        m_dim = m_index_flat->d;
        m_data_size = m_index_flat->ntotal;
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

    // @todo: To be removed after the experimental evaluation, since this function is not a common API
    std::pair<std::vector<std::vector<VidType>>, std::vector<std::vector<float>>> BatchSegmentKnnQuery(const std::vector<VectorDataType>& query_data_list, const size_t& query_k) {
        std::vector <std::vector<VidType>> answer_list;
        std::vector <std::vector<float>> answer_list_d;
        const size_t nq = query_data_list.size();
        const size_t answer_size = query_k * nq;

        faiss::idx_t *idx_list = new faiss::idx_t[answer_size];
        float *dist_list = new float[answer_size];
        float *xq = new float[m_dim * nq];

        for (size_t i = 0; i < nq; ++i) {
            for (size_t j = 0; j < m_dim; ++j) {
                xq[i * m_dim + j] = query_data_list[i].data[j];
            }
        }

        std::fill(idx_list, idx_list + answer_size, -1);
        std::fill(dist_list, dist_list + answer_size, 0.0);

        m_index_flat->search(nq, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i = 0; i < answer_size; ++i) {
                if (idx_list[i] >= 0 && dist_list[i] >= 0)
                    dist_list[i] = std::sqrt(dist_list[i]);
            }
        }

        for (size_t i = 0; i < nq; ++i) {
            std::vector <VidType> knn_list;
            std::vector<float> knn_list_d;
            for (size_t j = 0; j < query_k; ++j) {
                size_t idx = idx_list[i * query_k + j];
                knn_list.emplace_back(m_data[idx].vid);
                knn_list_d.emplace_back(dist_list[i * query_k + j]);
            }
            answer_list.emplace_back(knn_list);
            answer_list_d.emplace_back(knn_list_d);
        }

        delete[] idx_list;
        delete[] dist_list;
        delete[] xq;

        return std::make_pair(answer_list, answer_list_d);
    }

    size_t IndexSize() const override {
        size_t index_size = 0;

        // VectorDataType* m_data;
        // faiss::IndexFlat* m_index_flat;
        index_size += sizeof(void *) * 2;

        // size_t code_size;
        // std::vector<uint8_t> codes;
        index_size += sizeof(size_t);
        index_size += sizeof(uint8_t) * m_index_flat->codes.size() + sizeof(size_t);

        // std::vector<float> cached_l2norms;
        if (std::is_same<Distance, EuclideanDistance>::value ||
            std::is_same<Distance, EuclideanSquareDistance>::value) {
            faiss::IndexFlatL2 *index_tmp = dynamic_cast<faiss::IndexFlatL2 *>(m_index_flat);
            if (index_tmp)
                index_size += sizeof(float) * index_tmp->cached_l2norms.size() + sizeof(size_t);
        }

        return index_size;
    }

    size_t DataSize() const override {
        return m_data_size;
    }

    void SaveIndex(const std::string& index_path) const override {
        faiss::write_index(m_index_flat, index_path.c_str());
    }

private:
    void m_KnnQuery(const VectorDataType& query_data, const size_t& query_k, std::vector<VectorDataType>& answer_object, std::vector<VectorDimensionType>& answer_dist) {
        answer_object.clear();
        answer_dist.clear();

        faiss::idx_t *idx_list = new faiss::idx_t[query_k];
        float *dist_list = new float[query_k];
        float *xq = new float[m_dim];

        for (size_t j = 0; j < m_dim; ++j) {
            xq[j] = query_data.data[j];
        }

        std::fill(idx_list, idx_list + query_k, -1);
        std::fill(dist_list, dist_list + query_k, 0.0);

        m_index_flat->search(1, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i = 0; i < query_k; ++i) {
                if (idx_list[i] >= 0 && dist_list[i] >= 0)
                    dist_list[i] = std::sqrt(dist_list[i]);
            }
        }

        for (size_t i = 0; i < query_k; ++i) {
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

        faiss::idx_t *idx_list = new faiss::idx_t[answer_size];
        float *dist_list = new float[answer_size];
        float *xq = new float[m_dim * nq];

        for (size_t i = 0; i < nq; ++i) {
            for (size_t j = 0; j < m_dim; ++j) {
                xq[i * m_dim + j] = query_data_list[i].data[j];
            }
        }

        std::fill(idx_list, idx_list + answer_size, -1);
        std::fill(dist_list, dist_list + answer_size, 0.0);

        m_index_flat->search(nq, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i = 0; i < answer_size; ++i) {
                if (idx_list[i] >= 0 && dist_list[i] >= 0)
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

    void m_ClearIndexFlat() {
        if (m_index_flat != NULL)
            delete m_index_flat;
        m_index_flat = NULL;
        m_data = NULL;
    }

    faiss::IndexFlat *m_index_flat;
    const VectorDataType *m_data = NULL;
    int m_dim;
    size_t m_data_size;
};

#endif // INDEX_LINEARSCAN_INDEX_HPP