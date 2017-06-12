#include "mymatr.h"
#include "mytex.h"
#include <cstdio>
#include <jsoncpp/json/json.h> // or jsoncpp/json.h , or json/json.h etc.
int main () {
  myVector empty(0) ;
  size_t k=0, n=0, dep_mask=0, hash_rad=3, hash_max=20, itype=0 ;
  rand_t d1=-3, d2=3 ;
  char cdep_mask[32] ;
  
//  cin >> k >> n >> cdep_mask >> d1 >> d2 >> hash_rad >> hash_max >> itype ;

  if (!k) k = rand()%3+4 ;
  if (!n) n = 11 - k ;
      
  if (!dep_mask) {
	int mm = min(k, n-1) ;
	int dd2 = (1<<n) - 1 ;
	int dd1 = ( 1<<(n-1) ) | (1<<(n-3) ) ;

	if (!itype) itype = rand()%3 +1 ;
	if (itype==3) itype = (rand()%2) ? itype : 2 ;

	if (itype == 3) { // no solutions
  	  do {
		dep_mask = myrand (dd1, dd2) ;
	  } while (wt(dep_mask, n)<mm-2 || wt(dep_mask, n)>mm || !(dep_mask%2) ) ;
	}
	else if (itype == 1) { // one solution
	  dep_mask = ((1<<mm)-1) << ( n-mm ) ; 
	}
	else if (itype == 2) { //infinite solutions
  	  do {
		dep_mask = myrand (dd1, dd2) ;
	  } while (wt(dep_mask,n)<mm-2 || wt(dep_mask,n)>=mm || dep_mask%2) ;
	}
  }

	myMatrix A (k, n) ;
	A.gen_step(dep_mask, d1, d2) ;
	myMatrix B = A.hash(hash_rad, hash_max) ;
  Json::Value val ;
  val["name"]="Метод Гаусса" ;
  string statement, solution ;
  statement+="Исследовать систему методом Гаусса:\n" ;
	statement+=tex_equation_string(B) ;
	B.gaus1() ;
	
	myMatrix B_short = B ; //added
	B_short.del_col (B_short[0].size() - 1) ; // added
	pair<cnt_vect_t, string> temp=B_short.gaus1() ; //added
  statement+=temp.second ;
	statement+= "\n%" + string(itoa(k)) + ' ' + string(itoa(n));
  stringstream ss ;
	comment_matrix (B,"step",ss) ;
  statement+=ss.str() ;
	statement+= "%cols=\n%" ;
	
	solution+= "Ответ:\n\n" ;
	
	size_t rk_sys = B_short.main_vect.size() ; //added
	
	solution+= "ранг ситстемы: " + string(itoa(rk_sys)) + "\n\n" ; //changed
	solution+= "Главные неизвестные:\n" ;
	
	for (int i=0; i<rk_sys; i++) {
	  if (i) solution+= ", " ;
	  solution+="$x_" + string(itoa(B_short.main_vect[i]+1)) + "$" ;
	  statement+=string(itoa(B_short.main_vect[i]+1)) + ' ' ;
	}
	statement+= "\n%rk=\n%"  + string(itoa(rk_sys)) ;
	solution+= "\n\nСтупенчатый вид:\n$$\n" ;
	solution+=tex_matrix_string(B) ;
	solution+= "$$\n" ;
	//comment_matrix (A, "easy", cerr) ;
  val["solution"]=solution ;
  val["statement"]=statement ;
  cout << val << endl ;
  return 0 ;
}
