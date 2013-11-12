#ifndef OMP_PARALLEL_REDUCE
#define OMP_PARALLEL_REDUCE

#include "pfunc/pfunc.hpp"

namespace pfunc {

template <typename PFuncInstanceType,
          typename ReduceExecutable,
          typename SpaceType,
          typename PartitionerType>
struct omp_parallel_reduce : pfunc::virtual_functor {
  public:
  typedef typename PFuncInstanceType::taskmgr TaskMgrType;
  typedef typename PFuncInstanceType::task TaskType;

  private:
  /**
   * Used for applicatoring the function --- we need this so that we can comply
   * with PFunc's function type.
   */
  struct applicator : pfunc::virtual_functor {
    ReduceExecutable func;
    SpaceType space;
  
    applicator (const ReduceExecutable& func, const SpaceType& space) :
      func(func), space(space) {}
  
    void operator()() { func (space); }
  };

  SpaceType space;
  ReduceExecutable& func;
  TaskMgrType& taskmgr;
  
  public:
  omp_parallel_reduce (SpaceType space,
                       ReduceExecutable& func,
                       TaskMgrType& taskmgr) : 
    space(space), func(func), taskmgr(taskmgr) {}

  void operator()() {
    /* Figure out the number of threads */
    unsigned int size;
    pfunc::get_num_threads (taskmgr, size);
    
    /* Prepare all the tasks for launch */
    TaskType sub_tasks[size];
    std::auto_ptr<applicator> funcs[size];
    for (unsigned int rank=0; rank<size; ++rank) {
      SpaceType subspace = PartitionerType::create 
          (space.begin(), space.end(), rank, size);
      funcs[rank].reset(new applicator(func.split(rank), subspace));
    }
    
    /* Launch the tasks */
    for (unsigned int rank=0; rank<size; ++rank) {
      pfunc::spawn (taskmgr, sub_tasks[rank], *funcs[rank]);
    }
    
    /* wait for them */
    pfunc::wait_all (taskmgr, sub_tasks, sub_tasks+size);
    
    /* aggregate the results */
    for (unsigned int rank=0; rank<size; ++rank) {
      func.join(funcs[rank]->func);
    }
  } 
};

} // namespace pfunc

#endif // OMP_PARALLEL_REDUCE
