package algorithms.search.global;

/**
 * Goal is to make an implementation of quantum annealing global
 * search, that is,
 * quantum stochastic optimization, related to subject
 * adiabatic quantum computation.
 * 
 * First reading the paper "An introduction to quantum annealing"
 * by Falco and Tamascelli (2009?).
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
 *
   "Exponential Enhancement of the Efficiency of Quantum Annealing
      by Non-Stoquastic Hamiltonians"
      by Nishimori* and Takada, 2017
      
    ?
    * QUBO:
    * https://github.com/alex1770/QUBO-Chimera
      Quadratic unconstrained binary optimization (QUBO) is a 
      pattern matching technique, common in machine learning applications. 
      QUBO is an NP hard problem.
          QUBO problems may sometimes be well-suited to algorithms aided 
          by quantum annealing.[1]
    
    Chimera graph, a rectangularly linked grid of bipartite graphs used 
    as the hardware connectivity structure in quantum annealers
 
 * as an exercise, 
 * explore
 *     http://www.quantumplayground.net/#/home
 * 
 * might impl
 * Shor:
 *    integer factorization in a quantum computer.
 *    The efficiency of Shor's algorithm is due to the efficiency of the 
 *    quantum Fourier transform, and modular exponentiation by repeated 
 *    squarings.
 * 
 * Grovers:
 *     https://pdfs.semanticscholar.org/834b/579cde622a51a1f1459319a48452e5d9760d.pdf
 * 
 * https://github.com/dwavesystems/
 * 
 * looking at dwave code...
 * https://github.com/dwavesystems/dw_sa_chi/tree/master/src
 * chimera_annealer dwave_sa(fields, coupler_starts, coupler_ends, coupler_values, seed);
 
 
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
    
           refs 9 and 8, emphasize that they don't reproduce the computationally
              untractable details of a quantum mechanical system,
              but instead
              simulate the ground state process of the Hamiltonian Hν as
              nu approaches 0.
              hence an efficient exploration of solution space is possible.
    
              the ground state estimation (psi_nu) 
                 uses time evolution from the real portion of the
                 shrodinger equation
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
                        makes N(t) transitions towards nearest neighbors
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
                        psi_nu(y) for every neighbor y of current solution x, 
                        the exploration results in high probability of
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
    
    if pursue using mc, should look further into project 
    espresso of sissa.it.
    and TurboRVB
    
    ?
    https://github.com/hadsed/pathintegral-qmc

    https://github.com/alex1770/QUBO-Chimera

Look for recent articles by these authors:
Inoue J-I 2005 Quantum Annealing and Related Optimization Methods ed A Das and B K Chakrabarti
(Berlin:Springer) p 171

   a quick recent review:
    https://arxiv.org/pdf/1310.1339.pdf

   ---------------
   Notes and excerpts from
      "Optimization using quantum mechanics: quantum annealing through 
       adiabatic evolution" 
       by Santoro and Tosatti

   ... challenging hard optimization problems,
        such as the random Ising model, 
        the travelling salesman problem and 
        Boolean satisfiability problems. 
    The techniques used to implement quantum annealing are either 
    deterministic Schrodinger’s evolutions for the toy models, 
    or path-integral Monte Carlo and Green’s function Monte Carlo approaches, 
    for the hard optimization problems. 
    The crucial role played by disorder and the associated
    non-trivial Landau–Zener tunnelling phenomena is discussed and emphasized.

    The description as a toy model is because the Hilbert space of the
    problem is restricted to low N (N < 20 to 30 for spin 1/2 problems, e.g.).

    To examine the performance of different algorithms better, the authors
    look in detail at the energy structure and barriers around local minima,
    starting with simple models.

    First, one-dimensional potentials, starting from a
    double-well potential, the simplest form of barrier. 
    -- compared quantum adiabatic Schrodinger ¨
       evolution, both in real and in imaginary time, 
       and its classical deterministic counterpart, i.e.,
       Fokker–Planck evolution [28]. 
       (NOTE: Fokker-Planck evolution a.k.a. the Kolomogorv Forward eqn.
       ...time evolution of the PDE of the velocity of a particle under 
       the influence of drag forces and random forces, as in Brownian motion_.
    -- also studied the performance of different stochastic 
       approaches, both classical Monte Carlo and path-integral 
       Monte Carlo. 
    The authors found ~ equivalence of IT and RT Schrodinger annealing 
    justifying practical implementations of quantum annealing based on 
    imaginary-time quantum Monte Carlo schemes.

    the properties of the instantaneous spectrum of the Hamiltonian (26) 
    depend on the dimensionality D of the lattice, and on the 
    nearest-neighbour kinetic energy terms included 
    (adjacency matrix of the quantum walker)
   
    PIMC:
       Berne et al. combined clever thermal and quantum annealing schemes within a 
       path-integral framework to find the minimal energy configurations of
       protein models on the continuum 
       (with up to N = 46 monomers) [31, 32].
          the minimal energy configurations of Lennard-Jones clusters 
          of up to N = 100 atoms [33, 34]. 
       
       PIMC is summarized using 2 approaches:
         (1) disrete optimization problems, by using Ising case
         (2) continuum problems, using a particle in  potential

       simulates the equilibrium behavior of a system a finite temp T

       Edward–Anderson Ising glass in a transverse field:
           H = Hcl + Hkin = − summ_i_j (J_i_j * σ_i^2 + σ_j^2) - JMP(σ_i^x)

       this is a time-dependent Hamiltonian, and goal is to find min states
       as a function of time by turning off the transverse field (t).
       Can simulate the thermodynamics of a fixed positive T at a fixed JMP(t)
       by sampling the quantum partition function:

       The authors again quote the PIMC references above for interesting detailed
       solutions specific to their problems, 
       then present a simplified model to continue exploration:
          PIMC-QA scheme whereby T is untouched during a QA simulation 
          and the number of Trotter replicas P is also kept fixed, 
          in order to more clearly discriminate possibly genuine quantum 
          effects from thermal ones. 

          PIMC-QA applied to combinatorial optimization problems:
          
          define a combinatorial optimization as the task of 
          minimizing any given cost function that depends on 
          variables w/ discrete values.
          it is straightforward to map problems over the search 
          for the ground state of some Hamiltonian depending on Potts 
          (or Ising) spin degrees of freedom [7, 55]. 

          the number of possible configurations of a very small 
          32 × 32 square lattice Ising model is of order 10308, 
          while the number of electrons in the universe is of 
          the order 10^308

          The authors find better success for the PIMC-QA and GFMC models 
          than for the CA models.
          W.R.T. 2 MC-QA models, 
          the PIMC model difficulties are finite temperature T, 
          possible sampling problems for the action, and difficulties with 
          the Trotter break-up, so there may be alternative models that
          arise in PIMC specifically that are further improvements.
          the GFMC difficulty is that it is dependent upon
          a good guess of trial variational wavefunctions.

   ----------------
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
