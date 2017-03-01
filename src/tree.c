#include "tree.h"
#include "Lmer.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

extern long int alloc_mem;
extern long int nnodes;

Node* new_node(){
   int i;
   Node* newnode = (Node *) malloc(sizeof(Node)); 
   for (i = 0; i < ACGT; i++ ){
      newnode->children[i] = NULL;
   }
   return newnode;
}

static Node**  pool_ptr = NULL;
static int pool_count = 0;
static int pool_index = 0 ;

Node* get_new_pool(){
   Node * pool;
   if(pool_count % POOL_NUM == 0) {
      pool_ptr = realloc(pool_ptr,sizeof(Node*)*( pool_count + POOL_NUM));
      alloc_mem += sizeof(Node*)*(POOL_NUM);
   } 
   pool = malloc( sizeof(Node) * BUF_NODE );
   alloc_mem +=( sizeof(Node) * BUF_NODE );
   pool_ptr[pool_count++] = pool; 
   return pool;
}

Node* new_node_buf(){
   int i; 
   static Node * pool = 0 ;
   if (! pool_index ){
        pool = get_new_pool();
        pool_index = BUF_NODE ;
   }
   Node * newnode = pool ++ ;
   for (i = 0; i < ACGT; i++ ){
      newnode->children[i] = NULL;
   }
   nnodes++;
    
   pool_index-- ;
   return newnode;
}

void  free_tree(Node* tree){
   if (tree != NULL){
      for (int i = 0; i< ACGT; i++){
         if(tree->children[i] !=  NULL)
           free_tree( (tree->children[i]));
      }
   free(tree);
   }
}


void free_all_nodes(){
   int i;
   int N = (pool_count/POOL_NUM+1)*POOL_NUM;
   for ( i = 0; i < N ; i++){
      if(pool_ptr[i] != NULL){
          free(pool_ptr[i]);
          alloc_mem -= sizeof(Node) * BUF_NODE;
          nnodes -= BUF_NODE; 
      }
   }
   nnodes += pool_index; 
   free(pool_ptr);
   alloc_mem -= sizeof(Node *)*N; 
   pool_ptr = NULL;
   pool_count = 0; 
   pool_index = 0; 
}

void insert_Lmer(Node *tree, char *Lmer, int L){
  int i = 0; 
  Node* current = tree; // current and tree point to the same place
  for (i = 0; i < L ; i++ ){
      if ( (int) Lmer[i] > 3) break; // ignore N's 
      if (current->children[ (int) Lmer[i]] == NULL ){
         current->children[ (int) Lmer[i]] = new_node_buf();
      }
      current = current -> children[ (int)  Lmer[i]];
  }  
}

void construct_tree(Node *tree, char *sequence, int N, int L){
   int i;
   Lmer_sLmer(sequence,N);
   //Run over all arrays
   for ( i = 0; i < (N - L + 1); i++){
      insert_Lmer(tree, sequence+i,L);
   }
}

bool check_path(Node *tree, char *Lmer ,int L, int Lread){
  int i,j; 
  L = min(L, Lread);
  for (i = 0 ; i< Lread - L + 1; i++){
    Node* current = NULL;
    current = tree; 
    for (j = 0; j< L; j++){
        if (current->children[ (int) Lmer[i + j]] == NULL ){
           return false;
        }
        else{
           current = current -> children[ (int) Lmer[i + j]];
        }
    }
  }
  return true; 
}



