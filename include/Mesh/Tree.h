#ifndef ASN_TREE_H
#define ASN_TREE_H

#include "../Config.h"

namespace Asn
{
    template <typename DataType, Int NUM_RED_CHILD>
    class Tree;

    template <typename DataType, Int NUM_RED_CHILD>
    class Tree
    {
    public:
        struct Node
        {
            const DataType data;
            const Int level;

            Node *parent;
            List<Node *> child;

            explicit Node(const DataType &data, Node *parent = nullptr);
            Node(const Node &node) = delete;
            ~Node() = default;

            void set_child(const List<DataType> &data);
            void del_child();

            bool is_red() const;
            bool is_green() const;

            Node &operator=(const Node &node) = delete;
        };

    private:
        static void copy(Node *&dst, const Node *const src, Node *parent = nullptr);

    private:
        List<Node *> _root;

    public:
        const List<Node *> &root;

    public:
        Tree();
        Tree(const Tree &tree);
        ~Tree();

        Node *add_root_node(const DataType &data);

        void clear();

        Tree &operator=(const Tree &tree);
    };
};

#endif
