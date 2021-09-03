"""
Comparison and optimization of model spectra to data.
"""
import logging
logger = logging.getLogger('Inference')

import os,sys
import operator
import numpy
from numpy import logical_and, logical_not
import dadi
from dadi import Misc, Numerics
from scipy.special import gammaln
import scipy.optimize
from dadi.Inference import *
from dadi.Numerics import _cached_projection


#: Stores thetas
_theta_store = {}
#: Counts calls to object_func
_counter = 0
#: Returned when object_func is passed out-of-bounds params or gets a NaN ll.
_out_of_bounds_val = -1e8


def _object_func_marginals(params, data_vec, model_func, pts, 
                 lower_bound=None, upper_bound=None, 
                 verbose=0, multinom=True, flush_delay=0,
                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
                 output_stream=sys.stdout, store_thetas=False):
    """
    Objective function for optimization, when we optimize using multiple marginal SFSs.
    """
    
    nmarginals= len(data_vec)
    #print "nmarginals in marginal_optimization"
    #print nmarginals
    
    global _counter
    _counter += 1

    #if nmarginals <2:
    #	print "error: number of marginals less than two, but optimization function for multiple marginals is used!"
    # 	return dadi.Inference._out_of_bounds_val

    # Deal with fixed parameters
    params_up = dadi.Inference._project_params_up(params, fixed_params)

    # Check our parameter bounds
    if lower_bound is not None:
        for pval,bound in zip(params_up, lower_bound):
            if bound is not None and pval < bound:
                #print "failure in bounds!, pval<lower_bound"
                return -dadi.Inference._out_of_bounds_val/ll_scale
    if upper_bound is not None:
        for pval,bound in zip(params_up, upper_bound):
            if bound is not None and pval > bound:
                return -dadi.Inference._out_of_bounds_val/ll_scale
    
    
    all_ns = [data_vec[marg_num].sample_sizes for marg_num in range(nmarginals)]
    #print "in marginal_optimization, all_ns is"
    #print all_ns
    
         
    all_args = [params_up, all_ns] + list(func_args)
    # Pass the pts argument via keyword, but don't alter the passed-in 
    # func_kwargs
    func_kwargs = func_kwargs.copy()
    func_kwargs['pts'] = pts
    
    #print all_args
    #print func_kwargs
    all_sfs = model_func(*all_args, **func_kwargs)
    #this supposes that the two thetas are equal. This should be verified in the end! 
    if multinom:
	result=numpy.sum([ll_multinom(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])
    else:
        result = numpy.sum([ll(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])

    if store_thetas:
        global _theta_store
        dadi.Inference._theta_store[tuple(params)] = numpy.mean([optimal_sfs_scaling(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])
    
    # Bad result
    if numpy.isnan(result):
        result = dadi.Inference._out_of_bounds_val
    
    if (verbose > 0) and (_counter % verbose == 0):
        param_str = 'array([%s])' % (', '.join(['%- 12g'%v for v in params_up]))
        output_stream.write('%-8i, %-12g, %s%s' % (_counter, result, param_str,
                                                   os.linesep))
        Misc.delayed_flush(delay=flush_delay)

    return -result/ll_scale

warningIssued=False
def _object_func_marginals_coarse(coarsenings,log=False):
	def _object_func_marginals_c(params, data_vec, model_func, pts, 
	                 lower_bound=None, upper_bound=None, 
	                 verbose=0, multinom=True, flush_delay=0,
	                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
	                 output_stream=sys.stdout, store_thetas=False,nmarginals=2):
	    """
	    Objective function for optimization, when we optimize using multiple marginal SFSs.
	    """
	    #print "data vec is"
	    #print data_vec.shape
	    
	    
	    global _counter
	    _counter += 1
	
	    if nmarginals <2 and not warningiIssued:
	    	print "Warning: number of marginals less than two, but optimization function for multiple marginals is used!"
	    	warningIssued=True
	    	#return dadi.Inference._out_of_bounds_val
	
	    # Deal with fixed parameters
	    params_up = dadi.Inference._project_params_up(params, fixed_params)
	
	    # Check our parameter bounds
	    if lower_bound is not None:
	        for pval,bound in zip(params_up, lower_bound):
	            if bound is not None and pval < bound:
	                return -dadi.Inference._out_of_bounds_val/ll_scale
	    if upper_bound is not None:
	        for pval,bound in zip(params_up, upper_bound):
	            if bound is not None and pval > bound:
	                return -dadi.Inference._out_of_bounds_val/ll_scale
	    
	    
	    all_ns = [data_vec[marg_num].sample_sizes for marg_num in range(nmarginals)]
	    #print "in marginal_optimization, all_ns is"
	    #print all_ns
	    
	         
	    all_args = [params_up, all_ns] + list(func_args)
	    # Pass the pts argument via keyword, but don't alter the passed-in 
	    # func_kwargs
	    func_kwargs = func_kwargs.copy()
	    func_kwargs['pts'] = pts
	    all_sfs = model_func(*all_args, **func_kwargs)
	    allcoarse=[coarsen.flatten(coarsen.split(all_sfs[i],coarsenings[i])) for i in range(len(all_sfs))]
	    #this supposes that the two thetas are equal. This should be verified in the end! 
	    if multinom:
		result=numpy.sum([ll_multinom(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])
	    else:
	        result = numpy.sum([ll(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])
	
	    if store_thetas:
	        global _theta_store
	        dadi.Inference._theta_store[tuple(params)] = numpy.mean([optimal_sfs_scaling(all_sfs[marg_num], data_vec[marg_num]) for marg_num in range(nmarginals)])
	    
	    # Bad result
	    if numpy.isnan(result):
	        result = dadi.Inference._out_of_bounds_val
	    
	    if (verbose > 0) and (_counter % verbose == 0):
	        param_str = 'array([%s])' % (', '.join(['%- 12g'%v for v in params_up]))
	        output_stream.write('%-8i, %-12g, %s%s' % (_counter, result, param_str,
	                                                   os.linesep))
	        Misc.delayed_flush(delay=flush_delay)
	
	    return -result/ll_scale
	if not log:	
	    return _object_func_marginals_c
	else:	
	    def _object_func_marginals_c_log(log_params, *args, **kwargs):
	         """
	         Objective function for optimization in log(params).
	         """
	         return _object_func_marginals_c(numpy.exp(log_params), *args, **kwargs)
	return _object_func_marginals_c_log 
		  










def _object_func_marginals_log(log_params, *args, **kwargs):
    """
    Objective function for optimization in log(params).
    """
    return _object_func_marginals(numpy.exp(log_params), *args, **kwargs)
    
def _object_func_marginals_coarse_log(log_params, *args, **kwargs):
    """
    Objective function for optimization in log(params).
    """
    return _object_func_marginals_coarse(numpy.exp(log_params), *args, **kwargs)    
    
    
    
def optimize_log(p0, data, model_func, pts, lower_bound=None, upper_bound=None,
                 verbose=0, flush_delay=0.5, epsilon=1e-3, 
                 gtol=1e-5, multinom=True, maxiter=None, full_output=False,
                 func_args=[], func_kwargs={}, fixed_params=None, ll_scale=1,
                 output_file=None,nmarginals=1):
    """
    Optimize log(params) to fit model to data using the BFGS method.

    This optimization method works well when we start reasonably close to the
    optimum. It is best at burrowing down a single minimum.

    Because this works in log(params), it cannot explore values of params < 0.
    It should also perform better when parameters range over scales.

    p0: Initial parameters.
    data: Spectrum with data.
    model_function: Function to evaluate model spectrum. Should take arguments
                    (params, (n1,n2...), pts)
    lower_bound: Lower bound on parameter values. If not None, must be of same
                 length as p0.
    upper_bound: Upper bound on parameter values. If not None, must be of same
                 length as p0.
    verbose: If > 0, print optimization status every <verbose> steps.
    output_file: Stream verbose output into this filename. If None, stream to
                 standard out.
    flush_delay: Standard output will be flushed once every <flush_delay>
                 minutes. This is useful to avoid overloading I/O on clusters.
    epsilon: Step-size to use for finite-difference derivatives.
    gtol: Convergence criterion for optimization. For more info, 
          see help(scipy.optimize.fmin_bfgs)
    multinom: If True, do a multinomial fit where model is optimially scaled to
              data at each step. If False, assume theta is a parameter and do
              no scaling.
    maxiter: Maximum iterations to run for.
    full_output: If True, return full outputs as in described in 
                 help(scipy.optimize.fmin_bfgs)
    func_args: Additional arguments to model_func. It is assumed that 
               model_func's first argument is an array of parameters to
               optimize, that its second argument is an array of sample sizes
               for the sfs, and that its last argument is the list of grid
               points to use in evaluation.
               Using func_args.
               For example, you could define your model function as
               def func((p1,p2), ns, f1, f2, pts):
                   ....
               If you wanted to fix f1=0.1 and f2=0.2 in the optimization, you
               would pass func_args = [0.1,0.2] (and ignore the fixed_params 
               argument).
    func_kwargs: Additional keyword arguments to model_func.
    fixed_params: If not None, should be a list used to fix model parameters at
                  particular values. For example, if the model parameters
                  are (nu1,nu2,T,m), then fixed_params = [0.5,None,None,2]
                  will hold nu1=0.5 and m=2. The optimizer will only change 
                  T and m. Note that the bounds lists must include all
                  parameters. Optimization will fail if the fixed values
                  lie outside their bounds. A full-length p0 should be passed
                  in; values corresponding to fixed parameters are ignored.
                  For example, suppose your model function is 
                  def func((p1,f1,p2,f2), ns, pts):
                      ....
                  If you wanted to fix f1=0.1 and f2=0.2 in the optimization, 
                  you would pass fixed_params = [None,0.1,None,0.2] (and ignore
                  the func_args argument).
    ll_scale: The bfgs algorithm may fail if your initial log-likelihood is
              too large. (This appears to be a flaw in the scipy
              implementation.) To overcome this, pass ll_scale > 1, which will
              simply reduce the magnitude of the log-likelihood. Once in a
              region of reasonable likelihood, you'll probably want to
              re-optimize with ll_scale=1.
    """
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout
    #print "in opt,"
    #print data.shape
    args = (data, model_func, pts, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params, 
            ll_scale, output_stream)
    if nmarginals==1:
    	object_fun=dadi.Inference._object_func_log
    else:
    	object_fun=_object_func_marginals_log


    p0 = dadi.Inference._project_params_down(p0, fixed_params)
    outputs = scipy.optimize.fmin_bfgs(object_fun, 
                                       numpy.log(p0), epsilon=epsilon,
                                       args = args, gtol=gtol, 
                                       full_output=True,
                                       disp=False,
                                       maxiter=maxiter)
    xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag = outputs
    xopt = dadi.Inference._project_params_up(numpy.exp(xopt), fixed_params)

    if output_file:
        output_stream.close()

    if not full_output:
        return xopt
    else:
        return xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflag    


def optimize_log_fmin(p0, data, model_func, pts, 
                      lower_bound=None, upper_bound=None,
                      verbose=0, flush_delay=0.5, 
                      multinom=True, maxiter=None, 
                      full_output=False, func_args=[], 
                      func_kwargs={},
                      fixed_params=None, output_file=None,nmarginals=1):
    """
    Optimize log(params) to fit model to data using Nelder-Mead. 

    This optimization method make work better than BFGS when far from a
    minimum. It is much slower, but more robust, because it doesn't use
    gradient information.

    Because this works in log(params), it cannot explore values of params < 0.
    It should also perform better when parameters range over large scales.

    p0: Initial parameters.
    data: Spectrum with data.
    model_function: Function to evaluate model spectrum. Should take arguments
                    (params, (n1,n2...), pts)
    lower_bound: Lower bound on parameter values. If not None, must be of same
                 length as p0. A parameter can be declared unbound by assigning
                 a bound of None.
    upper_bound: Upper bound on parameter values. If not None, must be of same
                 length as p0. A parameter can be declared unbound by assigning
                 a bound of None.
    verbose: If True, print optimization status every <verbose> steps.
    output_file: Stream verbose output into this filename. If None, stream to
                 standard out.
    flush_delay: Standard output will be flushed once every <flush_delay>
                 minutes. This is useful to avoid overloading I/O on clusters.
    multinom: If True, do a multinomial fit where model is optimially scaled to
              data at each step. If False, assume theta is a parameter and do
              no scaling.
    maxiter: Maximum iterations to run for.
    full_output: If True, return full outputs as in described in 
                 help(scipy.optimize.fmin_bfgs)
    func_args: Additional arguments to model_func. It is assumed that 
               model_func's first argument is an array of parameters to
               optimize, that its second argument is an array of sample sizes
               for the sfs, and that its last argument is the list of grid
               points to use in evaluation.
    func_kwargs: Additional keyword arguments to model_func.
    fixed_params: If not None, should be a list used to fix model parameters at
                  particular values. For example, if the model parameters
                  are (nu1,nu2,T,m), then fixed_params = [0.5,None,None,2]
                  will hold nu1=0.5 and m=2. The optimizer will only change 
                  T and m. Note that the bounds lists must include all
                  parameters. Optimization will fail if the fixed values
                  lie outside their bounds. A full-length p0 should be passed
                  in; values corresponding to fixed parameters are ignored.
    (See help(dadi.Inference.optimize_log for examples of func_args and 
     fixed_params usage.)
    """
    #print p0	
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout
	
    args = (data, model_func, pts, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params, 1.0,
            output_stream)
    #if nmarginals==1:
    #	object_fun=dadi.Inference._object_func_log
    #else:
    object_fun=_object_func_marginals_log
    
    p0 = dadi.Inference._project_params_down(p0, fixed_params)
    #print "optimizing!"
    
    #print object_fun
    #print numpy.log(p0)
    #print object_fun(p0,data,model_func,pts, lower_bound=lower_bound,upper_bound=upper_bound,verbose=0,multinom=multinom,flush_delay=flush_delay,func_args=func_args,func_kwargs=func_kwargs,fixed_params=fixed_params, ll_scale=1,output_stream=sys.stdout)
              
    outputs = scipy.optimize.fmin(object_fun, numpy.log(p0), args = args,
                                  disp=False, maxiter=maxiter, full_output=True)
    xopt, fopt, iter, funcalls, warnflag = outputs
    xopt = dadi.Inference._project_params_up(numpy.exp(xopt), fixed_params)

    if output_file:
        output_stream.close()

    if not full_output:
        return xopt
    else:
        return xopt, fopt, iter, funcalls, warnflag 
    
    
def optimize_log_fmin_coarse(p0, coarsenings,data, model_func, pts, 
                      lower_bound=None, upper_bound=None,
                      verbose=0, flush_delay=0.5, 
                      multinom=True, maxiter=None, 
                      full_output=False, func_args=[], 
                      func_kwargs={},
                      fixed_params=None, output_file=None,nmarginals=1):
    """
    Optimize log(params) to fit model to data using Nelder-Mead. 

    This optimization method make work better than BFGS when far from a
    minimum. It is much slower, but more robust, because it doesn't use
    gradient information.

    Because this works in log(params), it cannot explore values of params < 0.
    It should also perform better when parameters range over large scales.

    p0: Initial parameters.
    data: Spectrum with data.
    model_function: Function to evaluate model spectrum. Should take arguments
                    (params, (n1,n2...), pts)
    lower_bound: Lower bound on parameter values. If not None, must be of same
                 length as p0. A parameter can be declared unbound by assigning
                 a bound of None.
    upper_bound: Upper bound on parameter values. If not None, must be of same
                 length as p0. A parameter can be declared unbound by assigning
                 a bound of None.
    verbose: If True, print optimization status every <verbose> steps.
    output_file: Stream verbose output into this filename. If None, stream to
                 standard out.
    flush_delay: Standard output will be flushed once every <flush_delay>
                 minutes. This is useful to avoid overloading I/O on clusters.
    multinom: If True, do a multinomial fit where model is optimially scaled to
              data at each step. If False, assume theta is a parameter and do
              no scaling.
    maxiter: Maximum iterations to run for.
    full_output: If True, return full outputs as in described in 
                 help(scipy.optimize.fmin_bfgs)
    func_args: Additional arguments to model_func. It is assumed that 
               model_func's first argument is an array of parameters to
               optimize, that its second argument is an array of sample sizes
               for the sfs, and that its last argument is the list of grid
               points to use in evaluation.
    func_kwargs: Additional keyword arguments to model_func.
    fixed_params: If not None, should be a list used to fix model parameters at
                  particular values. For example, if the model parameters
                  are (nu1,nu2,T,m), then fixed_params = [0.5,None,None,2]
                  will hold nu1=0.5 and m=2. The optimizer will only change 
                  T and m. Note that the bounds lists must include all
                  parameters. Optimization will fail if the fixed values
                  lie outside their bounds. A full-length p0 should be passed
                  in; values corresponding to fixed parameters are ignored.
    (See help(dadi.Inference.optimize_log for examples of func_args and 
     fixed_params usage.)
    """
    if output_file:
        output_stream = file(output_file, 'w')
    else:
        output_stream = sys.stdout

    args = (data, model_func, pts, lower_bound, upper_bound, verbose,
            multinom, flush_delay, func_args, func_kwargs, fixed_params, 1.0,
            output_stream)
    if nmarginals==1:
    	object_fun=_object_func_log
    else:
    	object_fun=_object_func_marginals_coarse_log(coarsenings)
    	
    p0 = dadi.Inference._project_params_down(p0, fixed_params)
    outputs = scipy.optimize.fmin(object_fun, numpy.log(p0), args = args,
                                  disp=False, maxiter=maxiter, full_output=True)
    xopt, fopt, iter, funcalls, warnflag = outputs
    xopt = dadi.Inference._project_params_up(numpy.exp(xopt), fixed_params)

    if output_file:
        output_stream.close()

    if not full_output:
        return xopt
    else:
        return xopt, fopt, iter, funcalls, warnflag     



def from_data_dict(data_dict, pop_ids, projections, mask_corners=True,
                       polarized=True,print_successrate=False):
        """
        Spectrum from a dictionary of polymorphisms.

        pop_ids: list of which populations to make fs for.
        projections: list of sample sizes to project down to for each
                     population.
        polarized: If True, the data are assumed to be correctly polarized by 
                   `outgroup_allele'. SNPs in which the 'outgroup_allele'
                   information is missing or '-' or not concordant with the
                   segregating alleles will be ignored.
                   If False, any 'outgroup_allele' info present is ignored,
                   and the returned spectrum is folded.

        The data dictionary should be organized as:
            {snp_id:{'segregating': ['A','T'],
                     'calls': {'YRI': (23,3),
                                'CEU': (7,3)
                                },
                     'outgroup_allele': 'T'
                    }
            }
        The 'calls' entry gives the successful calls in each population, in the
        order that the alleles are specified in 'segregating'.
        Non-diallelic polymorphisms are skipped.
        
        print_successrate: prints the proportion of misses sites due to insufficient alleles in any population, or lack of polarization. Such information is important for estimating the effective sequencing length.  
        """
        failed=0
        Npops = len(pop_ids)
        fs = numpy.zeros(numpy.asarray(projections)+1)
        for snp, snp_info in data_dict.items():
            # Skip SNPs that aren't triallelic.
            if len(snp_info['segregating']) != 2:
                continue

            allele1,allele2 = snp_info['segregating']
            if not polarized:
                # If we don't want to polarize, we can choose which allele is
                # derived arbitrarily since we'll fold anyways.
                outgroup_allele = allele1
            elif 'outgroup_allele' in snp_info\
               and snp_info['outgroup_allele'] != '-'\
               and snp_info['outgroup_allele'] in snp_info['segregating']:
                # Otherwise we need to check that it's a useful outgroup
                outgroup_allele = snp_info['outgroup_allele']
            else: 
                # If we're polarized and we didn't have good outgroup info, skip
                # this SNP.
                failed+=1
                continue
    
            # Extract the allele calls for each population.
            allele1_calls = numpy.asarray([snp_info['calls'][pop][0]
                                            for pop in pop_ids])
            allele2_calls = numpy.asarray([snp_info['calls'][pop][1]
                                            for pop in pop_ids])
            # How many chromosomes did we call successfully in each population?
            successful_calls = allele1_calls + allele2_calls
            
            # Which allele is derived (different from outgroup)?
            if allele1 == outgroup_allele:
                derived_calls = allele2_calls
            elif allele2 == outgroup_allele:
                derived_calls = allele1_calls
    
            # To handle arbitrary numbers of populations in the fs, we need
            # to do some tricky slicing.
            slices = [[numpy.newaxis] * len(pop_ids) for ii in range(Npops)]
            for ii in range(len(pop_ids)):
                slices[ii][ii] = slice(None,None,None)
        
            # Do the projection for this SNP.
            pop_contribs = []
            iter = zip(projections, successful_calls, derived_calls)
            for pop_ii, (p_to, p_from, hits) in enumerate(iter):
                contrib = _cached_projection(p_to,p_from,hits)[slices[pop_ii]]
                pop_contribs.append(contrib)
            addamount=reduce(operator.mul, pop_contribs)
            print addamount.sum()
            fs += addamount
        fsout = dadi.Spectrum.Spectrum(fs, mask_corners=mask_corners, 
                         pop_ids=pop_ids)
        if polarized:
            return fsout
        else:
            return fsout.fold()


    
