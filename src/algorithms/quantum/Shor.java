package algorithms.quantum;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.GreatestCommonDenominator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.Random;
import thirdparty.ardeleanasm.gates.EGateTypes;
import thirdparty.ardeleanasm.gates.ExpModGate;
import thirdparty.ardeleanasm.gates.GateProducer;
import thirdparty.ardeleanasm.gates.GatesAbstractFactory;
import thirdparty.ardeleanasm.gates.IGate;
import thirdparty.ardeleanasm.quantum.exception.RegisterOverflowException;
import thirdparty.ardeleanasm.quantum.utils.QRegisterOperations;
import thirdparty.ardeleanasm.qubits.QRegister;
import thirdparty.ardeleanasm.qubits.Qubit;

/**
 * NOT READY FOR USE YET
 * 
 Implementation of Shor's factoring algorithm to find the prime factors of 
 integer N.
 
 following pseudocode in wikipedia
 https://en.wikipedia.org/wiki/Shor%27s_algorithm
 
 Note that part of the code needs to use a simulated quantum computer on a 
 classical computer.
 The contrast in space complexity is quickly large:
    the state of an n-bit system requires n bits,
    but an n-qubit on a classical computer requires 
    pow(2, n) complex coefficients.
 
 In quantum computation, instead of the probable states of all bits 
 adding up to 1,
 allowed operations are unitary matrices, which are effectively rotations 
 and they preserve the sum of the squares equals one, 
 the Euclidean or L2 norm.
 Also, because the states are represented by spin, 
 quantum computations are reversible.
 
*/
public class Shor {
    
     /*
    from wikipedia
    
    Shor's algorithm consists of two parts:
    -- A reduction, which can be done on a classical computer, 
       of the factoring problem to the problem of order-finding.
    -- then A quantum algorithm to solve the order-finding problem.
    
       1) Pick a random number a < N.
       2) Compute gcd(a, N). This may be done using the Euclidean algorithm.
       3) If gcd(a, N) ≠ 1, then this number is a nontrivial factor of N, 
          so we are done.
       4) Otherwise, use the period-finding subroutine (below) to find r, 
          the period of the following function:
              f(x)=a^{x}mod N
              smallest pos integer r for which f(x+r) = f(x)
              f(x + r)=a^{x + r}mod N = a^{x}mod N
       5) If r is odd, go back to step 1.
       6) If a^{r/2} ==  −1 (mod N), go back to step 1.
       7) gcd(a^{r/2} + 1, N) and gcd(a^{r/2} - 1, N) 
          are both nontrivial factors of N. We are done.
       For example: N=15,a=7,r=4
          gcd(7^{2}\pm 1,15)== gcd(49\pm 1,15)} 
              where gcd(48,15)=3 and gcd(50,15)=5
      
    *) Quantum part: Period-finding subroutine
    
      see details in article or in libquantum for this part
    
    For exact QFFT, reading
        https://arxiv.org/pdf/quant-ph/0301093.pdf
    
    standard QFFT is for n qbits and has order 2^n
    (Coppersmith [5] (see also Shor [11])
    
    ** R. B. Griffiths and C. Niu, 
    Semiclassical Fourier Transform for Quantum
    Computation, Phys. Rev. Lett. 76 (1996) pp.3228-3231 
        (also quant-ph/9511007)
    */
    
    private final Random rng;
    
    private final int number;
    
    public Shor(int number) {
               
        if (number < 6) {
            throw new IllegalArgumentException("number must be > 5");
        }
        
        long seed = System.currentTimeMillis();
        rng = Misc.getSecureRandom();
        rng.setSeed(seed);
        
        System.out.println("SEED=" + seed);
        
        this.number = number;
    }

    public Shor(int number, long seed) {
        
        if (number < 6) {
            throw new IllegalArgumentException("number must be > 5");
        }
        
        rng = Misc.getSecureRandom();
        rng.setSeed(seed);
        
        System.out.println("SEED=" + seed);
        
        this.number = number;
    }

    /**
     * NOT READY FOR USE YET
     * 
     * @return 
     */
    public int[] run() {
        
        //NOTE: could make this class extend QuantumAlgorithms

        int i;
  
        // try random GCD first (no quantum computer emulation)
        int[] result = randomGCD(number);
        assert(result != null);
        
        //3) If gcd(a, N) ≠ 1, then this number is a nontrivial factor of N, 
        //   so we are done.
        
        // counting the number of factors that are not 1, and are positive 
        int a = 0;
        int nFactors = 0;
        for (int r : result) {
            if (r > 1) {
                nFactors++;
                if (r > a) {
                    a = r;
                }
            }
        }
        if (nFactors > 1) {
            // temporarily commenting out while impl rest of algorithm
    //        return result;
        } else if (nFactors == 0) {
            throw new IllegalArgumentException("either " + number + " is not "
                + " valid input or there is an error in the algorithm");
        }
  
        System.out.println("factorization so far=" + Arrays.toString(result));
        System.out.println("random factor=" + a + " of number=" + number);
        
        /*
        4) Otherwise, use the period-finding subroutine (below) to find r, 
          the period of the following function:
              f(x)=a^{x}mod N
              smallest pos integer r for which f(x+r) = f(x)
              f(x + r)=a^{x + r}mod N = a^{x}mod N
        */
        
        int maxQ = 2 * number * number - 1;
        int log2Nsq = (int)Math.ceil(Math.log(maxQ)/Math.log(2));
        
        /*
        The quantum circuits used for this algorithm are custom designed for 
        each choice of N and each choice of the random a used in 
              f(x) = a*x mod N. 
        Given N, find Q = 2^q such that N^2 .leq. Q .lt. 2N^2, 
        which implies Q/r>N. The input and output qubit registers need to hold 
        superpositions of values from 0 to Q − 1, and so have q qubits each. 
        Using what might appear to be twice as many qubits as necessary 
        guarantees that there are at least N different x which produce the 
        same f(x), even as the period r approaches N/2.
        
        The circuit diagram uses symbols [H] for hadamard operations,
        U*a*2^n for _______ operations, and QFT quantum FT.
        
                                            ______
        |0> -[H]----------------------O-----|     |----
        ...                          ...    |     | ...
        |0> -[H]------------O---------|-----| QFT |----
        |0> -[H]---O--------|---------|-----|_____|----
        |1> -------|--------|---------|-------
                   |        |         |
                [Ua2^0]  [Ua2^1] [Ua2^n-1]
        */
                
        // ------- find the order r of a modulo N ------
        
        /*
        period finding
        step P1) apply Hadamard gates to all qubits in the input register. 
                 since Q is between N^2 and 2*N^2, need 2*N^2 qubits.
                 Q^(-1/2) * summation( |x> )
                            x=0 to Q-1
        */
        
        GatesAbstractFactory factory = GateProducer.getGateFactory();
        
        // number^2 <= Q <= 2*number^2
        // Q = 2^q --> q = log_2(Q)
        System.out.println("n qubits=" + log2Nsq);
        
        /*
        notes from "Basic concepts in quantum computation" by 
           Ekert, Hayden, and Inamori
        
        We shall assume that information is stored in the registers in binary 
        form.  For example, 
            the number 6 is represented by a register in state 
            |1> ⊗ |1> ⊗ |0>. 
        In more compact notation: |a> stands for the 
           tensor product |a_(n−1)> ⊗ |a_(n−2)>. . . |a_1> ⊗ |a_0>, 
           where a_i ∈ {0, 1},
        and it represents a quantum register
        prepared with the value 
            a = 2^(0)a_0 + 2^(1)a_1 + . . . 2^(n−1)a_(n−1). 
        There are 2^n states of this kind, representing all binary strings 
        of length n or numbers from 0 to 2^(n−1),
        and they form a convenient computational basis. In the following 
        a ∈ {0, 1}^n
        (a is a binary string of length n) implies that 
        |a> belongs to the computational basis.
        
        Thus a quantum register of size three can store individual numbers 
        such as 3 or 7,
            |0> ⊗ |1> ⊗ |1> ≡ |011> ≡ |3>, 
            |1> ⊗ |1> ⊗ |1> ≡ |111> ≡ |7>, 
        
        In the same quantum register of size 3 qubits.
            the superposition (1/√2) * (|0> + |1>).
        is 
            (1/√2) * (|0> + |1>)  ⊗  (1/√2) * (|0> + |1>)  ⊗  (1/√2) * (|0> + |1>)
            
        which in binary is:
            |000> + |001> + |010> + |011> + |100> + |101> + |110> + |111>
        
        and in decimal notation:
            |0> + |1> + |2> + |3> + |4> + |5> + |6> + |7>
        
        the later is the language the algorithm might be using that
        needs to be transposed into the qubit operations.
        
        The later is the decimal notation for a superposition of all states
        for a 3 qubit system.
        
        ...
        Any attempt to measure the state
            α|0> + β|1>
        results in |0> with probability |α|2, and |1> with probability |β|2.
        
        */
       
        
        // initialize the registers to a superposition of Q states.
        // this can be done by applying hadamard to all qubits
        // OR can use the quantum fourier transform.
        
        QRegister qReg = new QRegister(log2Nsq);
        qReg.initialize();
        
        QRegisterOperations regOps = QRegisterOperations.getInstance();
		IGate hGate = factory.getGate(EGateTypes.E_HadamardGate);	

        // tensor product of all qubits:
        Qubit superposition = regOps.entangle(qReg); 
        
        int[] targetPosition = new int[log2Nsq];
        for (int ii = 1; ii < log2Nsq; ++ii) {
            targetPosition[ii] = ii;
        }
        
        //double qCo = Math.pow(maxQ, -0.5);
        superposition = hGate.applyGate(superposition, targetPosition, null, log2Nsq);	

        //qReg.change(0, superposition);
        
        System.out.println("reg=" + qReg);
        System.out.println("sp=" + superposition);
        
        assert(superposition.isValid());
        
        /*
        period finding
        step P2)  Construct f(x) as a quantum function and apply it to the 
          above state, to obtain 

              Q^(-1/2) * summation over Q bits ( |x, f(x) >

          This is still a superposition of Q states. 

          f(x) = a^x mod N   
        
        */
        
        IGate eGate = new ExpModGate();
        int[] conditions = new int[] {a, number};
        
        //NOTE: paused here
        
        superposition = eGate.applyGate(superposition, targetPosition, 
            conditions, log2Nsq);
        
        assert(superposition.isValid());

        //qReg.change(0, superposition);
        
        System.out.println("2 sp=" + superposition);
        
        
        /*
        period finding
        step P3)  quantum fourier transform
        
        can be performed efficiently 
        on a quantum computer, with a particular decomposition into 
        a product of simpler unitary matrices. Using a simple 
        decomposition, the discrete Fourier transform on 
        2^(n) amplitudes can be 
        implemented as a quantum circuit consisting of only 
        O(n^2) Hadamard gates and controlled phase shift gates, 
        where n is the number of qubits.[1] This can be compared with 
        the classical discrete Fourier transform, which takes 
        O(n*2^n) gates (where n is the number of bits), 
        which is exponentially more than O(n^2). However, the quantum 
        Fourier transform acts on a quantum state, whereas the 
        classical Fourier transform acts on a vector, so not every 
        task that uses the classical Fourier transform can take 
        advantage of this exponential speedup.

        Can be implemented with unitary gates
             controlled phase gate and hadamard gate
        
        libquantum:
            for(i=width-1; i>=0; i--) {
                for(j=width-1; j>i; j--) {
                    quantum_cond_phase(j, i, reg);
                }
                quantum_hadamard(i, reg);
    
          where quantum_cond_phase
             float phi = (pi / ((MAX_UNSIGNED) 1 << (control - target)));
             complex z = cos(phi) + IMAGINARY * sin(phi)
             for(i=0; i<reg->size; i++) {
                reg->node[i].amplitude *= z;
             quantum_decohere(reg)
        
        libquantum is distributed under the terms of the GNU General Public
        License, version 3 (GPLv3), which is located in the file COPYING.

        Send inquiries, comments, bug reports, suggestions, patches, etc. to:
        libquantum@libquantum.de

        See also the libquantum website:
        http://www.libquantum.de/
        */
        
        
        /*
       5) If r is odd, go back to step 1.
       6) If a^{r/2} ==  −1 (mod N), go back to step 1.
       7) gcd(a^{r/2} + 1, N) and gcd(a^{r/2} - 1, N) 
          are both nontrivial factors of N. We are done.
       For example: N=15,a=7,r=4
          gcd(7^{2}\pm 1,15)== gcd(49\pm 1,15)} 
              where gcd(48,15)=3 and gcd(50,15)=5
      
        */
      
    return null;
}

    public int[] randomGCD(int N) {
        
        if (N < 6) {
            throw new IllegalArgumentException("N must be > 5");
        }
        
        /*
        1) Pick a random number a < N.
        2) Compute gcd(a, N). This may be done using the Euclidean algorithm.
        */
        
        TIntSet tried = new TIntHashSet();
        
        int[] result = null;
        
        int nFactors = 0;
        while (nFactors < 2) {
        
            int a = rng.nextInt(N - 2) + 2;
            while (tried.contains(a)) {
                a = rng.nextInt(N - 2) + 2;
            }
            tried.add(a);
        
            result = GreatestCommonDenominator.extendedEuclid(a, N);
        
            for (int r : result) {
                if (r > 1) {
                    nFactors++;
                }
            }
        }
        
        return result;
    }

    private String createZeroes(int n) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; ++i) {
            sb.append("0");
        }
        return sb.toString();
    }
}
