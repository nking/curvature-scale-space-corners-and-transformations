package algorithms.search.global;

/**
 * Goal is to make an implementation of quantum annealing glocal
 * search, that is,
 * quantum stochastic optimization.
 * 
 * First reading the paper "An introduction to quantum annealing"
 * by Falco and Tamascelli (year?).
 * 
 * It contains a review of simulating quantum computing,
 * quantum annealing optimization and then an improvement
 * of the later with a diffusion term.
 * 
 * The paper also mentions modern variants of quantum annealing
 * and that may be what is implemented here at the end of the
 * reading.
 
 * might also read in detail
 *     "Quantum path minimization: An efficient method for global optimization"
 *     by Liu and Berne, 2003
 */

/*
   Notes and excerpts while reading
   -------------------
 * Quantum Annealing, or Quantum Stochastic Optimization, is a
   classical randomized algorithm
 * 
 *  In simulated annealing the space of admissible solutions to a given 
 * optimization problem is visited by a temperature dependent random walk.
 * The cost function defines the potential energy profile of the solution 
 * space and thermal fluctuations prevent the exploration from getting stuck
 * in a local minimum. 
 * An opportunely scheduled temperature lowering (annealing), then, stabilizes
   the walk around a, hopefully global, minimum of the potential profile.
   The idea of using quantum, instead of thermal, jumps to explore the solution
   space of a given optimization problem was proposed 
   in "A numerical implementation
   of Quantum Annealing" by Apolloni, Cesa-Bianchi, and  Falco. 
   and in "Stochastic Processes" by Albeverio et al.
.  It was suggested by the behaviour of the stochastic process 
   q_ν associated with the ground state (state of minimal energy) of a 
   Hamiltonian of the form:
                v^2     
         H_v = ---- * (2nd deriv w.r.t. x) + V(x)
                 2      
 
 * section 3, outline
    - creating classical heuristics
    - admissible solutions to an optimization problem is
       defined by a temperature dependent random walk
       for simulated annealing
       - cost function has
           - potential energy of solution space
           - thermal fluctuations (to escape local minima)
              - specialy timed temperature lowering (annealing)
                to regulate the walk after the escape
              - QUANTUM replaces the thermal jump with tunneling
    - desribing the quantum jump:
       q_nu is associated w/ ground state of Hamiltonian
             v^2     
       H_v = ---- * (2nd deriv w.r.t. x) + V(x)
              2 
         -> the potential V encodes the cost function to be minimized
    
       Given ground state psi_nu of H_nu of a fixed nu,
       the stochastic process
       q_nu is built by "the ground state transformation",
         which is not yet defined here.
        
       the ground state transformation results in
                   
       H_v = -nu * L_nu + E
                        v^2     
           where L_v = ---- * (2nd deriv w.r.t. x) + b_nu
                         2 

           where b_v = (1/2) * (1st deriv w.r.t. x) * ln(psi_vi(x)^2)

           as nu approaches 0, q_nu will behave like a Markov chain
              having disrete and stalbe state space.
    
           refs 9 and 8, emphasize that they donot reporduce the computationally
              untractable details of a quantum mechanical system,
              but instead
              simulate the ground state process of the Hamiltonian Hν as
              nu approaches 0.
              hence an efficient exploration of solution space is possible.
    
              the ground state estimation (psi_nu) 
                 is result of a fake time evolution
                 from the real portion of the
                 shrodinger equation, which is just the heat equation.
                  (exp(-t*H_nu))
              the ground state (psi_nu) is then described as
                  time approaches infinity for
    
                    (1/alpha_t) * exp(-t*H_nu/h) | psi_nu(0)  \ = psi_nu
                                                 |            /
    
                    where alpha_t = /phi_0 | psi(0) \ * exp(-t*E_0)
                                    \      |        /
            
               and eigenstates are phi_0(= psi_nu), phi_1,...
    
               and those eigenvalues E_0 < E_1 < ... E_N
    
            in the limit of t -> inf:
               (1/alpha_t) * exp(-t*H_nu/h) | psi_nu(0)  \ 
                                            |            /
               = psi_nu
    
               = <limit t->inf | phi_0> 
                  + (1/<phi_0|psi(0)>)*exp(-t*(E_n-E_0)/hbar) | phi_n>
                  * <phi_n | psi(0)>
    
          the authors then use the Feyman-Kac formula to estimate
          the evolved state exp(-t*H_nu) * psi
               given init state psi
               -- x is defined as expected value of integral of
                    the potential V (==cost) along the stochastic
                    trajectories eps(tau).
                      each eps(tau) starts at x and from time 0 to t
                        make N(t) transitions towards nearest neighbors
                        using a uniform probability.
                        N(t) is poisson process of intensity nu 
                           (NLK: presumably they mean intensity is
                            proportional to nu unless they are using
                            nu as intensity instead of frequency)
                ------> ** the details of the sampling and 
                            ground state (psi) estimation are in
                            Ref. 8.
                        when that is implemented, 
                        given an estimate of 
                        psi_nu(y) for every neighbor y of current solution
                        x, the exploration results in high probability of
                        finding the nearest-neighbour solution with  
                        max estimated value of psi_nu. 
                      Then the decision rule on which neighbour to accept
                      is based on tne ensemble eps(tau).
    
    When implemented in a working computer program, 
    the procedure described above requires a variety of 
    approximations. 
        First of all, equation (7) reproduces the 
        ground-state only in the limit t → ∞, 
        corresponding to sample paths of infinite length. 
        The expected value of the right hand side of (7), moreover, 
        will be estimated by means of a finite size sample. 
    Finally, the number of neighbours of a configuration may be too 
    large to allow the estimation of ψν on the whole neighbourhood 
    of a given solution. The accuracy of the approximation depends, 
    for example, on the actual length νt of the sample paths, 
    the number of paths n per neighbour 
    and the dimension |Neigh| of the subset of the set of neighbours 
    of a given solution to consider at each move. 
    Some “engineering” is also in order [9] 
    (see the pseudo-code of Quantum Stochastic Optimization reported here): 
    for example, we can call a local optimization procedure every tloc 
    quantum transitions. In addition, if the search looks to be stuck 
    for in a local minimum, we can force a jump to another local minimum.
    
    Ref 8: "Quantum stochastic optimization" 1989 by Apolloni et al
    Reg 9: "A numerical implementation of Quantum Annealing" 1990 by
           Apolloni et al.
       
    
    Consider variants:
    Imaginary Time Quantum Monte Carlo (ITQMC) 
        [26, 45, 46, 38]
    Ref 26: "Minimization of the potential energy surface of Lennard—Jones 
        clusters by quantum optimization"
        by Gregor and Car 2005
    Ref 45: "Optimization using quantum mechanics: quantum annealing through 
        adiabatic evolution" 2006 by Santoro and Tosatti
    Ref 46: same as 45, publ in 2008 in 
        J. Phys. A: Math. Theor., 41:209801, 2008
    Ref 38: "On product, generic and random generic quantum satisfiability."
         by Laumann et al. 2009
    
    if pursue an using mc, should look further into project 
    espresso of sissa.it.
    and TurboRVB
    
    ?
    https://github.com/hadsed/pathintegral-qmc
    
 * 
 * @author nichole
 */
public class QuantumAnnealing {
    
    /*
    pseudocode from "An introduction to quantum annealing"
        by Falco and Tamascelli
    
    NOTE: still reading the variants of QA and QSO so might
    not implement following the pseudocode here.
    
    Looking at the pseudocode, it seems as if the success of 
    finding the true global minimum depends upon the initial
    conditions, that is, the global min solution would have to be
    within a connected distance of tunelling and local minima from
    the start conditions...so the MC approach is probably better
    for some problems at exploring state space, but haven't read
    details yet...
    
    --------------------------------------------
    procedure 1 QA
    
    Input : initial condition init; control param nu;
        duration t_max; tunnel time t_drill; local opt time t_loc

    t = 0
    eps = init
    v_min = cost(eps)
    while t < t_max do
       j = 0
       repeat
          i = 0
          repeat
             eps = QuantumTransition(eps, nu, t_max)
             if cost(eps) < v_min then
                v_min = cost(c)
                i = 0; j = 0;
             else
                i = j + 1
             endif
          until i > t_loc
          epsilon = Local Optimization(eps)
          if cost(eps) < v_min then
             v_min = cost(eps)
             j = 0
          endif
       until j < t_drill
       draw a trajectory of length nu*t_max and jump there
       Local Optimization(c)
    
    end while

    -------------------------------
    Procedure 2 Quantum Transitions
    
    Input: initial condition eps; chain length nu*t; set of 
        neighbors to estimate Neigh
    
    for all neighbour k in Neigh do
       estimate the wave function psi_nu(k)
    end for
    best = select a neighbour in Neigh w/ probability prop to psi)nu
    
    return best

    --------------------------------
    Procedure 3 Local Optimization
    
    Input: initial condition eps
       return the best solution found by any steepest descent strategy
    */
  
    public QuantumAnnealing() {    
    }
    
}
