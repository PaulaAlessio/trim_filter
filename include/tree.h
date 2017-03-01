// tree.h 
#ifndef TREE_H
#define TREE_H
#include "defines.h"



typedef struct node{
   struct node *children[ACGT];
}Node;



Node *new_node();

Node* get_new_pool();

Node *new_node_buf();


void free_tree(Node* node);

void free_all_nodes();

void insert_Lmer(Node *tree, char *Lmer,int L);

void construct_tree(Node *tree,char *sequence,  int N, int L);

bool check_path(Node *tree, char *Lmer, int L, int Lread);



#endif
