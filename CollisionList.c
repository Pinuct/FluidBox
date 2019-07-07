#include <stdlib.h>
#include "CollisionList.h"

void removeItem(Collision_t *collision);

CollisionList_t *CollisionList_create(){
    CollisionList_t *head = NULL;
    head = malloc(sizeof(CollisionList_t));

    head->listSize = 0;
    head->next = NULL;
    return head;
}

void CollisionList_addElement(CollisionList_t *list, int i, int j, float dx, float dy, float r2) {
    Collision_t **current = &(list->next);
    while (*current != NULL) {
        current = &((*current)->next);
    }

    /* now we can add a new variable */
    *current = malloc(sizeof(Collision_t));
    (*current)->element_i = i;
    (*current)->element_j = j;
    (*current)->distance_x = dx;
    (*current)->distance_y = dy;
    (*current)->r2 = r2;
    (*current)->next = NULL;

    //aumentar el tamaño
    list->listSize++;
}

Collision_t *CollisionList_getElement(CollisionList_t *list, int position){
    Collision_t *current = list->next;

    // si la lsta esta vacia
    if (current == NULL) {
            return NULL;
    }

    int i;
    for (i = 0; i < position-1; i++) {
        if (current == NULL) {
            return NULL;
        }
        current = current->next;
    }

    return current;
}

void CollisionList_delete(CollisionList_t *list){
    removeItem(list->next);
    free(list);
}

void removeItem(Collision_t *collision){
    if(collision->next != NULL){
        removeItem(collision->next);
    }
    free(collision);
}