#ifndef ASN_TREE_HPP
#define ASN_TREE_HPP

#include "Tree.h"

namespace Asn
{
    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD>::Node::Node(const DataType &data, Tree<DataType, NUM_RED_CHILD>::Node *parent)
        : data(data), level(parent == nullptr ? 0 : parent->level + 1), parent(parent) {}

    template <typename DataType, Int NUM_RED_CHILD>
    inline void Tree<DataType, NUM_RED_CHILD>::Node::set_child(const List<DataType> &data)
    {
        assert(child.empty());
        for (Int i = 0; i < static_cast<Int>(data.size()); i++)
        {
            child.push_back(new Node(data[i], this));
        }
        return;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline void Tree<DataType, NUM_RED_CHILD>::Node::del_child()
    {
        for (Int i = 0; i < static_cast<Int>(child.size()); i++)
        {
            child[i]->del_child();
            delete child[i];
        }
        child.clear();
        return;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline bool Tree<DataType, NUM_RED_CHILD>::Node::is_red() const
    {
        return this->parent == nullptr ? true : this->parent->child.size() == NUM_RED_CHILD;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline bool Tree<DataType, NUM_RED_CHILD>::Node::is_green() const
    {
        return !is_red();
    }

    template <typename DataType, Int NUM_RED_CHILD>
    void Tree<DataType, NUM_RED_CHILD>::copy(Node *&dst, const Node *const src, Node *parent)
    {
        if (src == nullptr)
        {
            return;
        }
        dst = new Node(src->data, parent);
        dst->child.resize(src->child.size(), nullptr);
        for (Int i = 0; i < static_cast<Int>(dst->child.size()); i++)
        {
            copy(dst->child[i], src->child[i], dst);
        }
        return;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD>::Tree() : root(_root) {}

    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD>::Tree(const Tree<DataType, NUM_RED_CHILD> &tree) : Tree()
    {
        _root.resize(tree.root.size(), nullptr);
        #pragma omp parallel for
        for (Int i = 0; i < static_cast<Int>(_root.size()); i++)
        {
            copy(_root[i], tree.root[i]);
        }
        return;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD>::~Tree()
    {
        clear();
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD>::Node *Tree<DataType, NUM_RED_CHILD>::add_root_node(const DataType &data)
    {
        _root.push_back(new Node(data));
        return _root.back();
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline void Tree<DataType, NUM_RED_CHILD>::clear()
    {
        #pragma omp parallel for
        for (Int i = 0; i < static_cast<Int>(_root.size()); i++)
        {
            _root[i]->del_child();
            delete _root[i];
        }
        _root.clear();
        return;
    }

    template <typename DataType, Int NUM_RED_CHILD>
    inline Tree<DataType, NUM_RED_CHILD> &Tree<DataType, NUM_RED_CHILD>::operator=(const Tree<DataType, NUM_RED_CHILD> &tree)
    {
        if (&tree != this)
        {
            clear();
            _root.resize(tree.root.size(), nullptr);
            #pragma omp parallel for
            for (Int i = 0; i < static_cast<Int>(_root.size()); i++)
            {
                copy(_root[i], tree.root[i]);
            }
        }
        return (*this);
    }
};

#endif
