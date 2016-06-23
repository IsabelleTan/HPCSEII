//
// Created by Isabelle Tan on 04-05-16.
//


#ifndef EXERCISE3_QUADTREE_H
#define EXERCISE3_QUADTREE_H

// Set the value_type
typedef double value_type;

/*
 * A datastructure for nodes in the KD-Tree
 */
struct Node
{
    int level, morton_id;       // Nr of subdivisions in the tree for this node and morton index of square of this node.
    int child_id;               // Array index of the first of the four child nodes.
    int part_start, part_end;   // Array indices of first and last particle inside this node.
    value_type mass, xcom, ycom;     // The mass of this node (sum of particle mass) and x,y position of center of mass.
};


// Build a tree
void build(const value_type* const x, const value_type* const y, const value_type* const mass, const int N, const int k, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, Node* tree, int depth);

// A function that splits a parent node recursively as long as there are more than k particles inside
void split(Node* parent, Node* tree, int depth, uint32_t* index, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, int k, int *newNodeIndex);

// Subdivide a parent node
void assignParticles(Node* parent, Node* children, int depth, uint32_t* index);

// Compute the center of mass and the total mass in 4 children nodes.
void centerOfMass(Node* children, value_type* xsorted, value_type* ysorted, value_type* mass_sorted);

// A function that prints all the attributes of a Node
void printNode(Node node);

#endif //EXERCISE3_QUADTREE_H
