#ifndef loop

/* C-style loop (last not included) */
#define loop(I,FROM,TO) for ((I)=(FROM); (I)<(TO); (I)++)

/* Pascal etc. style loop (last included) */
#define loopto(I,FROM,TO) for ((I)=(FROM); (I)<=(TO); (I)++)

/* loop over linked list.  Example:
   struct mylist_s {
     struct mylist_s *next; // next item; NULL terminates the list
     char text[8]; } 
   *head,*l;
   looplist (l,head) printf("%s\n",l->text);
*/
#define looplist(PTR,HEAD) for ((PTR)=(HEAD); (PTR); (PTR)=(PTR)->next)

#endif
