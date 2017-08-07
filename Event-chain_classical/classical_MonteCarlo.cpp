#include<cmath>
#include<iostream>
#include<cstdio>
#include<vector>
#include<boost/random.hpp>

const int N = 200; // the number of particles
const int D = 3; // the dimension of the system
const double L = 5.0; // the system size


inline double two(double one);
inline double dist(double position1, double position2);
inline double dist2(double position1, double position2);
inline double irrdist2(std::vector<double> *i, std::vector<double> *j, int direction);
inline double potential(double position, double residist2, std::vector<double> *j, int direction, double beta); // (*i)[x] denotes x-coordinate of i-th particle
inline double inv_potential(double potential, double residist2, double beta, int which);
void output(std::vector<std::vector<double> > *particles,int index);


int main(){
  std::vector<std::vector<double> > particles(N,std::vector<double>(D));
  int direction = 0;
  int collision = 0;
  
  const int seed = 9739;
  boost::mt19937 eng(seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random(eng, boost::uniform_real<>());

  const int MCS = 4000000;
  const int THERMALIZATION = 90000;
  const double SWEEPS = 1.0;
  double sweep_interval = 0.0;
  double interval = 10.5;
  int measurement = 1;
  double temperature = 1.5;
  double beta = 1.0/temperature;
  
  int i,j;
  // initialize //
  for(i = 0; i < N; ++i){
    for(j = 0; j < D; ++j){
      particles[i][j] = L*random();
    }
  }
  output(&particles,0);
  // check //
  /*particles[0][0] = 0.0;
  particles[0][1] = 0.0;
  particles[0][2] = 0.0;
  particles[1][0] = 0.0;
  particles[1][1] = 1.0;
  particles[1][2] = 0.0;
  double checkirrdist2 = irrdist2(&particles[0],&particles[1],0);
  double checkpotential = potential(0.1,checkirrdist2,&particles[1],0,beta);
  std::cout << "##----check----##" << std::endl;
  std::cout << "#potential " << checkpotential << std::endl;
  std::cout << "#inv_potential " << inv_potential(checkpotential,checkirrdist2,beta,1) << std::endl;
  for(double xxx = -5; xxx < 5; xxx = xxx+0.001){
    std::cout << xxx << " " << potential(xxx,checkirrdist2,&particles[1],0,beta) << " " << inv_potential(xxx,checkirrdist2,beta,1) << " " << inv_potential(xxx,checkirrdist2,beta,-1) << " " << -inv_potential(xxx,checkirrdist2,beta,1) << " " << -inv_potential(xxx,checkirrdist2,beta,-1) << std::endl;
    }*/


  int monte_carlo_step;
  int moving = 0;
  int target = 0;
  int next_moving;
  double total_move = 0.0;
  for(monte_carlo_step = 0; monte_carlo_step < MCS; ++monte_carlo_step){
    double lift = L;
    for(target = 0; target < N; ++target){
      if(target != moving){
	double Estar = random();
	while( Estar == 0 ) Estar = random();
	Estar = -log(Estar);
	
	double residist2 = irrdist2(&particles[moving],&particles[target],direction);
	double dist_dir = dist(particles[moving][direction],particles[target][direction]);
	if(residist2 >= cbrt(2)){
	  double originalE;
	  if(dist_dir < 0) originalE = potential(particles[target][direction], residist2, &particles[target], direction,beta);
	  else originalE = potential(particles[moving][direction], residist2, &particles[target], direction,beta);
	  Estar += originalE;
	  double tmp = L;
	  if(Estar < 0) tmp = -dist_dir +inv_potential(Estar, residist2, beta,1);
	  if(tmp < lift) {
	    lift = tmp;
	    next_moving = target;
	    if(tmp < 0) std::cout << "#CHECK01" << std::endl;
	  }
	}
	else{
	  double originalE;
	  double minimum = sqrt(cbrt(2)-residist2);
	  double tmp = L;
	  if(dist_dir < -minimum || (dist_dir > 0 && dist_dir < minimum)) originalE = potential(particles[target][direction]+minimum,residist2,&particles[target],direction,beta);
	  else originalE = potential(particles[moving][direction],residist2,&particles[target],direction,beta);
	  Estar += originalE;
	  if(dist_dir < 0){
	    if(Estar < potential(particles[target][direction],residist2,&particles[target],direction,beta)){
	      tmp = -dist_dir - inv_potential(Estar, residist2,beta,1);
	      if(tmp < 0){
		std::cout << "#CHECK02 : tmp = " << tmp <<  std::endl;
	      }
	    }else{
	      Estar -= potential(particles[target][direction],residist2,&particles[target],direction,beta);
	      if(Estar < 0){
		tmp = -dist_dir + inv_potential(Estar, residist2, beta,-1);
		if(tmp < 0) std::cout << "#CHECK03" << std::endl;
	      }
	    }
	  }else{
	    if(Estar < 0){
	      tmp = -dist_dir + inv_potential(Estar, residist2,beta,-1);
	      if(tmp < 0) std::cout << "#CHECK04" << std::endl;
	    }
	  }
	  if(tmp < lift){
	    lift = tmp;
	    next_moving = target;
	  }
	}
      }
    }
    if(sweep_interval > SWEEPS && measurement < 5000){
      output(&particles,measurement);
      sweep_interval -= SWEEPS;
      measurement++;
    }
    if(lift == L){
      int arbitrary = floor(random()*N);
      next_moving = arbitrary;
      lift = random()*L;
    }else collision++;
    total_move += lift;
    sweep_interval += lift;
    if(total_move > interval){
      particles[moving][direction] += total_move - interval;
      while(particles[moving][direction] > L) particles[moving][direction] -=L;
      total_move = 0.0;
      direction = (direction+1)%D;
      interval = random() * 100 + 4.0;
    }else{
      particles[moving][direction] += lift;
      while(particles[moving][direction] > L) particles[moving][direction] -=L;
    }
    moving = next_moving;
  }
  std::cout << "collision" << collision << std::endl;
  
  return 0;
}

inline double two(double one){
  return one*one;
}

inline double dist(double position1, double position2){
  double dist = position1 - position2;
  while(dist > L/2){
    dist -= L;
  }
  while(dist < -L/2){
    dist += L;
  }
  return dist;
}

inline double dist2(double position1, double position2){
  double dist = position1 - position2;
  while (dist > L/2){
    dist -= L;
  }
  while (dist < -L/2){
    dist += L;
  }
  return dist*dist;
}

inline double irrdist2(std::vector<double> *i, std::vector<double> *j, int direction){
  if(D==1) return 0;
  else{
    double tmp = 0.0;
    int iterator;
    for(iterator = 1; iterator < D; ++iterator) tmp += dist2((*i)[(direction+iterator)%D],(*j)[(direction+iterator)%D]);
    return tmp;
  }
}
								     
inline double potential(double position, double residist2, std::vector<double> *j, int direction, double beta){
  double tmp = dist2(position, (*j)[direction]) + residist2;
  tmp = tmp*tmp*tmp;
  tmp = 1.0/tmp;
  return beta*(two(tmp) - tmp);
}

inline double inv_potential(double potential, double residist2, double beta, int which){
  return sqrt(cbrt(2.0/(1.0+which*sqrt(1+4*potential/beta))) - residist2);
}

void output(std::vector<std::vector<double> > *particles,int index){
  FILE *fp;
  char filename[80];
  int i,j;
  sprintf(filename, "./data/data%09d.dat",index);
  fp = fopen(filename, "w");
  for(i = 0; i < N; ++i){
    for(j = 0; j < D; ++j){
      fprintf(fp, "%lf ", (*particles)[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}
