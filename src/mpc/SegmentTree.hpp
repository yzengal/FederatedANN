#ifndef MPC_SEGMENT_TREE_HPP
#define MPC_SEGMENT_TREE_HPP

#include <vector>
#include <iomanip>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <limits>

class SegmentTree {
public:
    SegmentTree(int m): m_silo_num(m) {
        Init(m);
    }

    void Init(int m) {
        m_silo_num = m;
        if (m_silo_num == 0) {
            m_pos_list.clear();
            m_id_list.clear();
            return ;
        }

        std::vector<int> vec(m_silo_num);
        for (int i=0; i<m_silo_num; ++i) {
            vec[i] = i;
        }
        std::random_shuffle(vec.begin(), vec.end());

        m_pos_list.resize(m_silo_num);
        for (int i=0; i<m_silo_num; ++i) {
            m_pos_list[vec[i]] = i;
        }

        m_id_list.clear();
        m_id_list.emplace_back(vec);
        for (int i=0; ; ++i) {
            std::vector<int> pre_level = m_id_list[i];
            size_t pre_level_size = pre_level.size();
            if (pre_level_size == 1) break;

            std::vector<int> next_level;
            int min_id;
            for (int j=0; j<pre_level_size; j+=2) {
                if (j+1 >= pre_level_size) {
                    min_id = pre_level[j];
                } else {
                    min_id = std::min(pre_level[j], pre_level[j+1]);
                }
                next_level.emplace_back(min_id);
            }
            m_id_list.emplace_back(next_level);
        }
    }

    std::vector<int> LevelAt(int level_id) const {
        if (level_id >= m_id_list.size()) {
            std::string error_message;
            error_message = std::string("Segment tree has ") + std::to_string(m_id_list.size()) 
                            + std::string(" levels only, where level_id is ") + std::to_string(level_id);
            throw std::invalid_argument(error_message);
        }
        return m_id_list[level_id];
    }

    void UpdateLevelAt(int level_id, const std::vector<int>& vec) {
        if (level_id >= m_id_list.size()) {
            std::string error_message;
            error_message = std::string("Segment tree has ") + std::to_string(m_id_list.size()) 
                            + std::string(" levels only, where level_id is ") + std::to_string(level_id);
            throw std::invalid_argument(error_message);
        }
        if (m_id_list[level_id].size() != vec.size()) {
            std::string error_message;
            error_message = std::string("In Segment tree, Level #") + std::to_string(level_id) + std::string(" should have ") 
                            + std::to_string(m_id_list[level_id].size()) + std::string(" instead of ") 
                            + std::to_string(vec.size()) + std::string(" elements");
            throw std::invalid_argument(error_message);
        }
        m_id_list[level_id] = vec;
    }

    void UpdateNodeAt(const int level_id, const int node_idx, const int key) {
        if (level_id >= m_id_list.size()) {
            std::string error_message;
            error_message = std::string("Segment tree has ") + std::to_string(m_id_list.size()) 
                            + std::string(" levels only, where level_id is ") + std::to_string(level_id);
            throw std::invalid_argument(error_message);
        }
        if (node_idx >= m_id_list[level_id].size()) {
            std::string error_message;
            error_message = std::string("In Segment tree, Level #") + std::to_string(level_id) + std::string(" should have ") 
                            + std::to_string(m_id_list[level_id].size()) + std::string(" elements, while node_idx is ") + std::to_string(node_idx);
            throw std::invalid_argument(error_message);
        }
        m_id_list[level_id][node_idx] = key;
    }

    //        0
    //    0       2
    //  0   1   2
    // 0 1 2 3 4 5
    void GetParentNodeCompare(const int cur_level_id, const int cur_node_idx, const int cur_silo_id, int& level_id, int& node_idx, int& cmp_silo_id) {
        if (cur_silo_id >= m_silo_num) {
            std::string error_message;
            error_message = std::string("Segment tree has ") + std::to_string(m_silo_num) 
                            + std::string(" leaf nodes, where silo_id is ") + std::to_string(cur_silo_id);
            throw std::invalid_argument(error_message);
        }
        if (cur_level_id+1 == this->Level()) {
           level_id = node_idx = cmp_silo_id = -1;
           return ;
        }

        level_id = cur_level_id + 1;
        node_idx = cur_node_idx >> 1;

        int cmp_node_idx = ((cur_node_idx%2)>0) ? (cur_node_idx-1) : (cur_node_idx+1);
        if (cmp_node_idx >= m_id_list[cur_level_id].size()) {
            cmp_silo_id = cur_silo_id;
        } else {
            cmp_silo_id = m_id_list[cur_level_id][cmp_node_idx];
        }
    }

    int LeafNodeIdx(const int silo_id) const {
        return m_pos_list[silo_id];
    }

    int Root() const {
        if (m_id_list.empty()) return -1;
        int level = m_id_list.size();
        return *m_id_list[level-1].rbegin();
    }

    int LeafLevel() const {
        return 0;
    }

    int Level() const {
        return m_id_list.size();
    }

    std::string to_string() const {
        std::stringstream ss;

        for (int i=m_id_list.size()-1; i>=0; --i) {
            ss << "Level #" << i << ":";
            for (auto id : m_id_list[i]) {
                ss << " " << id;
            }
            ss << "\n";
        }

        return ss.str();
    }

private:
    size_t m_silo_num;
    std::vector<std::vector<int>> m_id_list;
    std::vector<int> m_pos_list;
};

#endif // MPC_SEGMENT_TREE_HPP