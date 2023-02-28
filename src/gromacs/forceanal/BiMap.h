/*
    BiMap.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/07/28
    Description: Bidirection Map. DO NOT require Boost library.
*/

#ifndef SRC_GROMACS_FORCEANAL_BIMAP_H_
#define SRC_GROMACS_FORCEANAL_BIMAP_H_

// C++ STL
#include <map>

namespace ForceAnal
{
    template<typename T1, typename T2>
    class BiMap
    {
    public:
        BiMap()
        {}

        void add(T1& key, T2& value)
        {
            forwardMap.insert(std::make_pair<T1, T2>(key, value));
            backwardMap.insert(std::make_pair<T2, T1>(value, key));
        }

        void add_asym(T1& left, T2& middle, T1& right)
        {
            forwardMap.insert(std::make_pair<T1, T2>(left, middle));
            backwardMap.insert(std::make_pair<T2, T1>(middle, right));
        }

        void remove_left(T1 key_left)
        {
            auto it_left = forwardMap.erase(key_left);
            auto it_right = backwardMap.erase(it_left->second);
        }

        T2& get_left(T1& key_left)
        {
            return forwardMap.at(key_left);
        }

        const T2& get_left(const T1& key_left)
        {
            return forwardMap.at(key_left);
        }

        T1& get_right(T2& key_right)
        {
            return backwardMap.at(key_right);
        }

        const T1& get_right(const T2& key_right)
        {
            return backwardMap.at(key_right);
        }

    private:
        std::map<T1, T2> forwardMap;
        std::map<T2, T1> backwardMap;
    };
} // ForceAnal

#endif /* SRC_GROMACS_FORCEANAL_BIMAP_H_ */
