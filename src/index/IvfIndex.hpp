#ifndef INDEX_IVF_INDEX_HPP
#define INDEX_IVF_INDEX_HPP

#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <faiss/IndexFlat.h>
#include <faiss/IndexIVFFlat.h>
#include <faiss/MetricType.h>
#include <faiss/invlists/DirectMap.h>
#include <faiss/invlists/InvertedLists.h>
#include <faiss/index_io.h>
#include <typeinfo>

#include "BaseIndex.hpp"
#include "../utils/MetricType.hpp"
#include "../utils/BenchLogger.hpp"

template<typename Distance = EuclideanDistance>
class IvfIndex final : public BaseIndex<Distance> {
public:
    IvfIndex() : m_ivf(NULL), m_data(NULL), m_dim(0), m_data_size(0) {

    }

    IvfIndex(const std::vector <VectorDataType> &data_list) : IvfIndex() {
        Construct(data_list);
    }

    IvfIndex(const std::string &index_path, const std::vector <VectorDataType> &data_list) : IvfIndex() {
        Load(index_path, data_list);
    }

    ~IvfIndex() override {
        delete m_ivf;
    }

    void Construct(const std::vector <VectorDataType> &data_list) override {
        m_ClearIvf();
        if (data_list.empty()) return;

        m_data = data_list.data();
        m_dim = data_list[0].Dimension();
        m_data_size = data_list.size();

        faiss::IndexFlatL2 quantizer(m_dim);
        m_ivf = new faiss::IndexIVFFlat(&quantizer, m_dim, m_ivf_nlist, faiss::MetricType::METRIC_L2);

        float *xb = new float[m_dim * m_data_size];

        for (size_t i = 0; i < m_data_size; ++i) {
            for (size_t j = 0; j < m_dim; ++j) {
                xb[i * m_dim + j] = data_list[i].data[j];
            }
        }

        m_ivf->train(m_data_size, xb);
        m_ivf->add(m_data_size, xb);

        delete[] xb;
    }

    void Load(const std::string &index_path, const std::vector <VectorDataType> &data_list) override {
        faiss::Index *loadedIndex = faiss::read_index(index_path.c_str());
        m_ivf = dynamic_cast<faiss::IndexIVFFlat *>(loadedIndex);
        m_ivf->nprobe = 20;
        m_data = data_list.data();
        m_dim = m_ivf->d;
        m_data_size = m_ivf->ntotal;
        std::cout << "Load hnsw index, dim = " << m_dim << ", data_size = " << m_data_size << std::endl;
    }

    size_t IndexSize() const override {
        size_t index_size = 0;
        // faiss::IndexIVFFlat *m_ivf;
        // const VectorDataType *m_data;
        index_size += sizeof(void *) * 2;

        // int d;        ///< vector dimension
        // idx_t ntotal; ///< total nb of indexed vectors
        // bool verbose; ///< verbosity level
        // bool is_trained;
        // MetricType metric_type;
        // float metric_arg;
        index_size += sizeof(int);
        index_size += sizeof(faiss::idx_t);
        index_size += sizeof(bool) * 2;
        index_size += sizeof(faiss::MetricType);
        index_size += sizeof(float);

        //// IndexIVF
        // InvertedLists* invlists = nullptr;
        index_size += sizeof(size_t) * 2;
        index_size += sizeof(bool);
        index_size += sizeof(size_t);
        for (int i = 0; i < m_ivf_nlist; ++i) {
            index_size += m_ivf->invlists->list_size(i) * sizeof(faiss::idx_t) + sizeof(size_t);
            index_size += m_ivf->invlists->list_size(i) * sizeof(unsigned char) + sizeof(size_t);
        }
        // bool own_invlists = false;
        // size_t code_size = 0;
        // int parallel_mode = 0;
        // const int PARALLEL_MODE_NO_HEAP_INIT = 1024;
        // bool by_residual = true;
        index_size += sizeof(bool);
        index_size += sizeof(size_t);
        index_size += sizeof(int) * 2;
        index_size += sizeof(bool);
        // DirectMap direct_map;
        index_size += sizeof(faiss::DirectMap::Type);
        index_size += m_ivf->direct_map.array.size() * sizeof(faiss::idx_t) + sizeof(size_t);
        index_size += m_ivf->direct_map.hashtable.size() * sizeof(faiss::idx_t) * 2 + sizeof(size_t);
        // size_t nprobe = 1;    ///< number of probes at query time
        // size_t max_codes = 0; ///< max nb of codes to visit to do a query
        index_size += sizeof(size_t) * 2;
        // Index* quantizer = nullptr;
        index_size += sizeof(size_t);
        //index_size += sizeof(uint8_t) * m_ivf_quantizer->codes.size() + sizeof(size_t);
        // if (std::is_same<Distance, EuclideanDistance>::value ||
        //     std::is_same<Distance, EuclideanSquareDistance>::value) {
        //     faiss::IndexFlatL2 *index_tmp = dynamic_cast<faiss::IndexFlatL2 *>(m_ivf_quantizer);
        //     if (index_tmp)
        //         index_size += sizeof(float) * index_tmp->cached_l2norms.size() + sizeof(size_t);
        // }
        // size_t nlist = 0;
        // char quantizer_trains_alone = 0;
        // bool own_fields = false; ///< whether object owns the quantizer
        // Index* clustering_index = nullptr; **
        index_size += sizeof(void *) * 2;
        index_size += sizeof(size_t);
        index_size += sizeof(char);
        index_size += sizeof(bool);
        // ClusteringParameters cp; ///< to override default clustering params
        index_size += sizeof(int) * 2;
        index_size += sizeof(bool) * 5;
        index_size += sizeof(int) * 3;
        index_size += sizeof(size_t);
        return index_size;
    }

    size_t DataSize() const override {
        return m_data_size;
    }

    std::vector <VectorDataType> KnnQuery(const VectorDataType &query_data, const size_t &query_k) override {
        std::vector <VectorDataType> answer_object;
        std::vector <VectorDimensionType> answer_dist;

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
        faiss::write_index(m_ivf, index_path.c_str());
    }

    void aaaa() {
        std::cerr << "hhhh" << std::endl;
        std::string store_index = "../faissindex/sift_IVF.faissindex";
        faiss::write_index(this->m_ivf, store_index.c_str());
        std::cout << store_index << std::endl;
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

        m_ivf->search(1, xq, query_k, dist_list, idx_list);

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

        m_ivf->search(nq, xq, query_k, dist_list, idx_list);

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

    void m_ClearIvf() {
        if (m_ivf != NULL)
            delete m_ivf;
        m_ivf = NULL;
        m_data = NULL;
    }

    faiss::IndexIVFFlat *m_ivf;
    const VectorDataType *m_data;
    int m_dim;
    size_t m_data_size;
    const int m_ivf_nlist = 100;
};


#endif // INDEX_IVF_INDEX_HPP