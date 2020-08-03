#include "pcross.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define  epsilon  0.001

#define  P_LEFT       -1
#define  P_RIGHT       1
#define  P_ON          0
#define  P_DIF_SIDE   -1

#define PItem  struct PITEM

PItem {
   real   x;
   real   y;
   int    pos;
   real   dist;
   PItem* next;
};


PItem*
MakePItem( real x, real y )
{
   PItem* res = (PItem*)malloc( sizeof(PItem) );
   assert( res != NULL );

   res->x = x;
   res->y = y;
   res->pos = P_ON;
   res->dist = 0;
   res->next = NULL;

   return res;
}


void
ClearList( PItem* list )
{
   PItem* tmp;
   if (list == NULL)  return;

   tmp = list->next;
   list->next = NULL;
   while  (tmp != NULL) {
      list = tmp;
      tmp = tmp->next;
      free(list);
   }
}


int
crosspolygons_( int* pn1, real* x1, real* y1,
               int* pn2, real* x2, real* y2,
               real* resX, real* resY )
{
   PItem*  list;
   int     N, i;
   PItem*  tmp;
   PItem*  tmp1;
   int     count;
   real    normx, normy, rtmp;
   int     n1, n2;

   n1 = *pn1;
   n2 = *pn2;
   list = MakePItem(x2[0], y2[0]);
   tmp = list;
   N = n2;
   for ( count = 1 ; count < N ; count++ ) {
      tmp->next = MakePItem( x2[count], y2[count] );
      tmp = tmp->next;
   }
   tmp->next = list;

   tmp = list;
   for ( count = 0 ; count < n1 ; count++ ) {
      if  (N <= 2)  break;

      normx = y1[(count+1)%n1] - y1[count];
      normy = x1[count] - x1[(count+1)%n1];
      rtmp = sqrt( normx*normx + normy*normy );
      normx /= rtmp;
      normy /= rtmp;

      i = N;
      do {
         tmp->dist = ( tmp->x - x1[count] ) * normx +
                     ( tmp->y - y1[count] ) * normy;

         tmp->pos = (tmp->dist < -epsilon) ? P_LEFT :
                    (tmp->dist > epsilon ) ? P_RIGHT : P_ON;

         tmp = tmp->next;
      } while (--i);

      i = N;
      do {
         if (tmp->pos * tmp->next->pos == P_DIF_SIDE)  {
            rtmp = tmp->dist / tmp->next->dist;
            tmp1 = MakePItem( (tmp->x - rtmp * tmp->next->x)/(1-rtmp),
                              (tmp->y - rtmp * tmp->next->y)/(1-rtmp) );
            tmp1->next = tmp->next;
            N++;
            tmp->next = tmp1;
            tmp = tmp1->next;
           } else
            tmp = tmp->next;
      } while (--i);

      i = N;
      do {
         if (tmp->next->pos == P_RIGHT)  {
            tmp1 = tmp->next;
            tmp->next = tmp1->next;
            free( tmp1 );
            if (--N == 2 ) break;
           } else
            tmp = tmp->next;
      } while (--i);    
      list = tmp;
   }

   if (N <= 2) {
      ClearList( list );
      return 0;
   }

   for ( count = 0 ; count < N ; count++ ) {
      resX[count] = tmp->x;
      resY[count] = tmp->y;
      tmp = tmp->next;
   }

   ClearList( list );
   return N;
}
#ifndef __PCROSS_H
#define __PCROSS_H

typedef double real;

/* 
  Функция CrossPolygons вычисляет пересечение двух выпуклых многоугольников;
      n1(n2) -- кол-во вершин в первом(втором) многоугольнике;
      x1(x2) -- значения кординат x вершин первого(второго) многоугольника;
      y1(y2) -- значения кординат y вершин первого(второго) многоугольника;
   Возвращает:
      если пересечение пусто -- 0;
      иначе -- кол-во вершин в полученном многоугольнике.
                                 * * *
   Вершины должны следовать в таком порядке, что при обходе от меньших
     индексов к большим внутренность многоугольника была слева.
                                 * * *
   Массивы resX и resY должны быть выделены до вызова функции; 
     и дочтаточны для размещения результата.                             */
 
int crosspolygons_( int* pn1, real* x1, real* y1,
                   int* pn2, real* x2, real* y2,
                   real* resX, real* resY ); 

#endif
