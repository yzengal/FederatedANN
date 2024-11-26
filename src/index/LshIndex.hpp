#ifndef INDEX_LSH_INDEX_HPP
#define INDEX_LSH_INDEX_HPP

#include <iostream>
#include <vector>
#include <stdexcept>
#include <faiss/IndexLSH.h>
#include <faiss/index_io.h>
#include <faiss/VectorTransform.h>

#include "BaseIndex.hpp"
#include "../utils/MetricType.hpp"
#include "../utils/BenchLogger.hpp"

template<typename Distance = EuclideanDistance>
class LshIndex final : public BaseIndex<Distance> {
public:
    LshIndex() : m_lsh(NULL), m_data(NULL), m_dim(0), m_data_size(0) {

    }

    LshIndex(const std::vector<VectorDataType>& data_list) : LshIndex() {
        Construct(data_list);
    }

    LshIndex(const std::string& index_path, const std::vector<VectorDataType>& data_list) : LshIndex() {
        Load(index_path, data_list);
    }

    ~LshIndex() {
        delete m_lsh;
    }

    void Construct(const std::vector<VectorDataType>& data_list) override {
        m_ClearLsh();
        if (data_list.empty()) return;

        m_data = data_list.data();
        m_dim = data_list[0].Dimension();
        m_data_size = data_list.size();

        m_lsh = new faiss::IndexLSH(m_dim, m_lsh_nbits);

        float *xb = new float[m_dim * m_data_size];

        for (size_t i = 0; i < m_data_size; ++i) {
            for (size_t j = 0; j < m_dim; ++j) {
                xb[i * m_dim + j] = data_list[i].data[j];
            }
        }

        m_lsh->add(m_data_size, xb);
        delete[] xb;
    }

    void Load(const std::string& index_path, const std::vector<VectorDataType>& data_list) override {
        faiss::Index* loadedIndex = faiss::read_index(index_path.c_str());
        m_lsh = dynamic_cast<faiss::IndexLSH*>(loadedIndex);
        m_data = data_list.data();
        m_dim = m_lsh->d;
        m_data_size = m_lsh->ntotal;
        std::cout << "Load hnsw index, dim = " << m_dim << ", data_size = " << m_data_size << std::endl;
    }

    size_t IndexSize() const override {
        size_t index_size = 0;
        // faiss::IndexLSH *m_lsh;
        // const VectorDataType *m_data;
        index_size += sizeof(void *) * 2;
        // size_t code_size;
        // std::vector<uint8_t> codes;
        index_size += sizeof(size_t);
        index_size += sizeof(uint8_t) * m_lsh->codes.size() + sizeof(size_t);
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
        // int nbits;             ///< nb of bits per vector
        // bool rotate_data;      ///< whether to apply a random rotation to input
        // bool train_thresholds; ///< whether we train thresholds or use 0
        // std::vector<float> thresholds; ///< thresholds to compare with
        index_size += sizeof(int);
        index_size += sizeof(bool) * 2;
        index_size += m_lsh->thresholds.size() * sizeof(float) + +sizeof(size_t);
        // RandomRotationMatrix rrot; ///< optional random rotation
        index_size += sizeof(bool) * 2;
        index_size += m_lsh->rrot.A.size() * sizeof(float) + sizeof(size_t);
        index_size += m_lsh->rrot.b.size() * sizeof(float) + sizeof(size_t);
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

    void KnnQuery(const VectorDataType &query_data, const size_t &query_k, std::vector <VectorDataType> &answer_object,
                  std::vector <VectorDimensionType> &answer_dist) override {
        m_KnnQuery(query_data, query_k, answer_object, answer_dist);
    }

    std::vector <std::vector<VectorDataType>>
    BatchKnnQuery(const std::vector <VectorDataType> &query_data_list, const size_t &query_k) override {
        std::vector <std::vector<VectorDataType>> answer_object_list;
        std::vector <std::vector<VectorDimensionType>> answer_dist_list;

        m_BatchKnnQuery(query_data_list, query_k, answer_object_list, answer_dist_list);

        return answer_object_list;
    }

    void SaveIndex(const std::string &index_path) const override {
        std::cout << "Saving hnsw index..." << std::endl;
        faiss::write_index(m_lsh, index_path.c_str());
    }

private:
    void
    m_KnnQuery(const VectorDataType &query_data, const size_t &query_k, std::vector <VectorDataType> &answer_object,
               std::vector <VectorDimensionType> &answer_dist) {
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

#ifdef LOCAL_DEBUG
        PrintLine(__LINE__);
        std::cout << "k = " << query_k << ", query_data = " << query_data.to_string() << std::endl;
        std::vector<float> tmp_float_list(xq, xq+m_dim);
        PrintLine(__LINE__);
        std::cout << "m_dim = " << m_dim << ", xq = ";
        PrintVector<float>(tmp_float_list);
#endif

        m_lsh->search(1, xq, query_k, dist_list, idx_list);

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

    void m_BatchKnnQuery(const std::vector <VectorDataType> &query_data_list, const size_t &query_k,
                         std::vector <std::vector<VectorDataType>> &answer_object_list,
                         std::vector <std::vector<VectorDimensionType>> &answer_dist_list) {
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

        m_lsh->search(nq, xq, query_k, dist_list, idx_list);

        if (std::is_same<Distance, EuclideanDistance>::value) {
            for (size_t i = 0; i < answer_size; ++i) {
                if (idx_list[i] >= 0 && dist_list[i] >= 0)
                    dist_list[i] = std::sqrt(dist_list[i]);
            }
        }

        for (size_t i = 0; i < nq; ++i) {
            std::vector <VectorDataType> knn_list;
            std::vector <VectorDimensionType> knnd_list;
            for (size_t j = 0; j < query_k; ++j) {
                size_t idx = idx_list[i * query_k + j];
                knn_list.emplace_back(m_data[idx]);
                knnd_list.emplace_back(dist_list[i * query_k + j]);
            }
            answer_object_list.emplace_back(knn_list);
            answer_dist_list.emplace_back(knnd_list);
        }

        delete[] idx_list;
        delete[] dist_list;
        delete[] xq;
    }

    void m_ClearLsh() {
        if (m_lsh != NULL)
            delete m_lsh;
        m_lsh = NULL;
        m_data = NULL;
    }


    faiss::IndexLSH *m_lsh;
    const VectorDataType *m_data;
    int m_dim;
    size_t m_data_size;
    int m_lsh_nbits = 16;
};

#endif // INDEX_LSH_INDEX_HPP