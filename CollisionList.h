#include <stdio.h>

typedef struct collision{
    int element_i;
    int element_j;
    float distance_x;
    float distance_y;
    float r2;
    struct collision *next;
} Collision_t;

typedef struct collisionList{
    int listSize;
    Collision_t *next;
} CollisionList_t;

CollisionList_t *CollisionList_create();
void CollisionList_addElement(CollisionList_t *list, int i, int j, float dx, float dy, float r2);
Collision_t *CollisionList_getElement(CollisionList_t *list, int position);
void CollisionList_delete(CollisionList_t *list);