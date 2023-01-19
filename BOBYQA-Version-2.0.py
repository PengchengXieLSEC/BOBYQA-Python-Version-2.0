#BOBYQA-Python-Version-2.0 Copyright: Pengcheng Xie Connect: xpc@lsec.cc.ac.cn

from math import sqrt, trunc, atan, cos, sin


def bobyqa(N, NPT, X, XL, XU, RHOBEG, RHO, IPRINT, MAXFUN):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION X(*),XL(*),XU(*),W(*)
    #
    #     This subroutine seeks the least value of a function of many variables,
    #     by applying a trust region method that forms quadratic models by
    #     interpolation. There is usually some freedom in the interpolation
    #     conditions, which is taken up by minimizing the Frobenius norm of
    #     the change to the second derivative of the model, beginning with the
    #     zero matrix. The values of the variables are constrained by upper and
    #     lower bounds. The arguments of the subroutine are as follows.
    #
    #     N must be set to the number of variables and must be at least two.
    #     NPT is the number of interpolation conditions. Its value must be in
    #       the interval [N+2,(N+1)(N+2)/2 - 1]. Choices that exceed 2*N+1 are not
    #       recommed.
    #     Initial values of the variables must be set in X[1 - 1],X[2 - 1],...,X[N - 1]. They
    #       will be changed to the values that give the least calculated F.
    #     For I=1,2,...,N, XL[I - 1] and XU[I - 1] must provide the lower and upper
    #       bounds, respectively, on X[I - 1]. The construction of quadratic models
    #       requires XL[I - 1] to be strictly less than XU[I - 1] for each I. Further,
    #       the contribution to a model from changes to the I-th variable is
    #       damaged severely by rounding errors if XU[I - 1]-XL[I - 1] is too small.
    #     RHOBEG and RHO must be set to the initial and final values of a trust
    #       \region radius, so both must be positive with RHO no greater than
    #       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
    #       expected change to a variable, while RHO should indicate the
    #       accuracy that is required in the final values of the variables. An
    #       error return occurs if any of the differences XU[I - 1]-XL[I - 1], I=1,...,N,
    #       is less than 2*RHOBEG.
    #     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
    #       amount of printing. Specifically, there is no output if IPRINT=0 and
    #       there is output only at the return if IPRINT=1. Otherwise, each new
    #       value of RHO is printed, with the best vector of variables so far and
    #       the corresponding value of the objective function. Further, each new
    #       value of F with its variables are output if IPRINT=3.
    #     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
    #     The array W will be used for working space. Its length must be at least
    #       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
    #
    #     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
    #     F to the value of the objective function for the current values of the
    #     variables X[1 - 1],X[2 - 1],...,X[N - 1], which are generated automatically in a
    #     way that satisfies the bounds given in XL and XU.
    #
    #     Return if the value of NPT is unacceptable.
    #
    NP = N + 1
    if (NPT < N + 2 or NPT > ((N + 2) * NP) / 2):
        print('Return from BOBYQA because NPT is not in  the required interval')
        return
    #
    #     Partition the working space array, so that different parts of it can
    #     be treated separately during the calculation of BOBYQB. The partition
    #     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
    #     space that is taken by the last array in the argument list of BOBYQB.
    #
    SU = [0 for _ in range(N)]
    SL = [0 for _ in range(N)]
    #
    #     Return if there is insufficient space between the bounds. Modify the
    #     initial X if necessary in order to avoid conflicts between the bounds
    #     and the construction of the first quadratic model. The lower and upper
    #     bounds on moves from the updated X are set now, in the ISL and ISU
    #     partitions of W, in order to provide useful and exact information about
    #     components of X that become within distance RHOBEG from their bounds.
    #
    ZERO = 0.0e0
    for J in range(1, N+1):
        TEMP = XU[J - 1] - XL[J - 1]
        if (TEMP < RHOBEG + RHOBEG):
            print('Return from BOBYQA because one of the' +
                  ' differences XU[I - 1]-XL[I - 1] is less than 2*RHOBEG.')
            return
        SL[J - 1] = XL[J - 1] - X[J - 1]
        SU[J - 1] = XU[J - 1] - X[J - 1]
        if (SL[J - 1] >= -RHOBEG):
            if (SL[J - 1] >= ZERO):
                X[J - 1] = XL[J - 1]
                SL[J - 1] = ZERO
                SU[J - 1] = TEMP
            else:
                X[J - 1] = XL[J - 1] + RHOBEG
                SL[J - 1] = -RHOBEG
                SU[J - 1] = max(XU[J - 1] - X[J - 1], RHOBEG)
        elif (SU[J - 1] <= RHOBEG):
            if (SU[J - 1] <= ZERO):
                X[J - 1] = XU[J - 1]
                SL[J - 1] = -TEMP
                SU[J - 1] = ZERO
            else:
                X[J - 1] = XU[J - 1] - RHOBEG
                SL[J - 1] = min(XL[J - 1] - X[J - 1], -RHOBEG)
                SU[J - 1] = RHOBEG
    #
    #     Make the call of BOBYQB.
    #
    bobyqb(N, NPT, X, XL, XU, RHOBEG, RHO, IPRINT, MAXFUN, SL, SU)


def bobyqb(N, NPT, X, XL, XU, RHOBEG, RHO, IPRINT, MAXFUN, SL, SU):
    #          DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #      1     XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),
    #      2     SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
    #
    #     The arguments N, NPT, X, XL, XU, RHOBEG, RHO, IPRINT and MAXFUN
    #       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
    #     XBASE holds a shift of origin that should reduce the contributions
    #       from rounding errors to values of the model and Lagrange functions.
    #     XPT is a two-dimensional array that holds the coordinates of the
    #       interpolation points relative to XBASE.
    #     FVAL holds the values of F at the interpolation points.
    #     XOPT is set to the printlacement from XBASE of the trust region centre.
    #     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
    #     HQ holds the explicit second derivatives of the quadratic model.
    #     PQ contains the parameters of the implicit second derivatives of the
    #       quadratic model.
    #     BMAT holds the last N columns of H.
    #     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
    #       this factorization being ZMAT times ZMAT**T, which provides both the
    #       correct rank and positive semi-definiteness.
    #     NDIM is the first dimension of BMAT and has the value NPT+N.
    #     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
    #       All the components of every XOPT are going to satisfy the bounds
    #       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
    #       XOPT is on a constraint boundary.
    #     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
    #       vector of variables for the next call of CALFUN. XNEW also satisfies
    #       the SL and SU constraints in the way that has just been mentioned.
    #     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
    #       in order to increase the denominator in the updating of UPDATE.
    #     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
    #     VLAG contains the values of the Lagrange functions at a new point X.
    #       They are part of a product that requires VLAG to be of length NDIM.
    #     W is a one-dimensional array thatNP is used for working space. Its length
    #       must be at least 3*NDIM = 3*(NPT+N).
    #
    #     Set some constants.
    #
    HALF = 0.5e0
    ONE = 1.0e0
    TEN = 10.0e0
    TENTH = 0.1e0
    TWO = 2.0e0
    ZERO = 0.0e0
    CAUCHY = ZERO
    NP = N + 1
    NDIM = NPT + N
    NPTM = NPT - NP
    NH = (N * NP) / 2
    #
    # Init worspace
    #
    XBASE = [0 for _ in range(N)]
    XPT = [[0 for _ in range(N)] for _ in range(NPT)]
    FVAL = [0 for _ in range(NPT)]
    XOPT = [0 for _ in range(N)]
    GOPT = [0 for _ in range(N)]
    HQ = [0 for _ in range(trunc((N * NP) / 2))]
    PQ = [0 for _ in range(NPT)]
    BMAT = [[0 for _ in range(N)] for _ in range(NDIM)]
    ZMAT = [[0 for _ in range(NPT - NP)] for _ in range(NPT)]
    XNEW = [0 for _ in range(N)]
    XALT = [0 for _ in range(N)]
    D = [0 for _ in range(N)]
    VLAG = [0 for _ in range(NDIM)]
    W = [0 for _ in range(3 * NDIM)]
    #
    #     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    #     BMAT and ZMAT for the first iteration, with the corresponding values of
    #     of NF and KOPT, which are the number of calls of CALFUN so far and the
    #     index of the interpolation point at the trust region centre. Then the
    #     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
    #     less than NPT. GOPT will be updated if KOPT is different from KBASE.
    #
    (X, XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NF, KOPT) \
        = prelim(N, NPT, X, XL, XU, RHOBEG, IPRINT, MAXFUN, XBASE, XPT,
                 FVAL, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU)
    XOPTSQ = ZERO
    for I in range(1, N+1):
        XOPT[I - 1] = XPT[KOPT - 1][I - 1]
        XOPTSQ = XOPTSQ + XOPT[I - 1] ** 2
    FSAVE = FVAL[0]
    if (NF < NPT):

        if (IPRINT > 0):
            print(['Return from BOBYQA because CALFUN has been' +
                  ' called MAXFUN times.'])
        if (FVAL[KOPT - 1] <= FSAVE):
            for I in range(1, N+1):
                X[I - 1] = min(max(XL[I - 1], XBASE[I - 1] +
                               XOPT[I - 1]), XU[I - 1])
                if (XOPT[I - 1] == SL[I - 1]):
                    X[I - 1] = XL[I - 1]
                if (XOPT[I - 1] == SU[I - 1]):
                    X[I - 1] = XU[I - 1]

            F = FVAL[KOPT - 1]

        if (IPRINT >= 1):
            print(
                f'At the return from BOBYQA. Number of function values = {NF}')

            print(f'Least value of F = {F}')
            print(f'The corresponding X is: {X[: N]}')
        return

    KBASE = 1
    #

    #     Complete the settings that are required for the iterative procedure.
    #
    RHO = RHOBEG
    DELTA = RHO
    NRESC = NF
    NTRITS = 0
    DIFFA = ZERO
    DIFFB = ZERO
    ITEST = 0
    NFSAV = NF
    #
    #     Update GOPT if necessary before the first iteration and after each
    #     call of RESCUE that makes a call of CALFUN.
    #
    flag = 20
    while (flag):
        if (flag == 20):
            if (KOPT != KBASE):
                IH = 0
                for J in range(1, N+1):
                    for I in range(1, J+1):
                        IH = IH + 1
                        if (I < J):
                            GOPT[J - 1] = GOPT[J - 1] + \
                                HQ[IH - 1] * XOPT[I - 1]
                        GOPT[I - 1] = GOPT[I - 1] + HQ[IH - 1] * XOPT[J - 1]

                if (NF > NPT):
                    for K in range(1, NPT+1):
                        TEMP = ZERO
                        for J in range(1, N+1):
                            TEMP = TEMP + XPT[K - 1][J - 1] * XOPT[J - 1]
                        TEMP = PQ[K - 1] * TEMP
                        for I in range(1, N+1):
                            GOPT[I - 1] = GOPT[I - 1] + \
                                TEMP * XPT[K - 1][I - 1]
            flag = 60
            #

            #     Generate the next point in the trust region that provides a small value
            #     of the quadratic model subject to the constraints on the variables.
            #     The integer NTRITS is set to the number "trust region" iterations that
            #     have occurred since the last "alternative" iteration. If the length
            #     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
            #     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
            #
        elif (flag == 60):
            (XNEW, D, W[:N], W[N:2 * N], W[2*N:3 * N], W[3 * N:4 * N],
             W[4 * N:5 * N], DSQ, CRVMIN) = \
                trsbox(N, NPT, XPT, XOPT, GOPT, HQ,
                       PQ, SL, SU, DELTA, XNEW, D, W[:N], W[N:2 * N],
                       W[2*N:3 * N], W[3 * N:4 * N], W[4 * N:5 * N])
            DNORM = min(DELTA, sqrt(DSQ))
            if (DNORM < HALF * RHO):
                while (True):
                    NTRITS = -1
                    DISTSQ = (TEN * RHO) ** 2
                    if (NF <= NFSAV + 2):
                        flag = 650
                        break

                    #
                    #     The following choice between labels 650 and 680 deps on whether or
                    #     not our work with the current RHO seems to be complete. Either RHO is
                    #     decreased or termination occurs if the errors in the quadratic model at
                    #     the last three interpolation points compare favourably with predictions
                    #     of likely improvements to the model within distance HALF*RHO of XOPT.
                    #
                    ERRBIG = max([DIFFA, DIFFB, DIFFC])
                    FRHOSQ = 0.125e0 * RHO * RHO
                    if (CRVMIN > ZERO and ERRBIG > FRHOSQ * CRVMIN):
                        flag = 650
                        break

                    BDTOL = ERRBIG / RHO
                    for J in range(1, N+1):
                        BDTEST = BDTOL
                        if (XNEW[J - 1] == SL[J - 1]):
                            BDTEST = W[J - 1]

                        if (XNEW[J - 1] == SU[J - 1]):
                            BDTEST = -W[J - 1]

                        if (BDTEST < BDTOL):
                            CURV = HQ((J + J * J) / 2)
                            for K in range(1, NPT+1):
                                CURV = CURV + PQ[K - 1] * \
                                    XPT[K - 1][J - 1] ** 2

                            BDTEST = BDTEST + HALF * CURV * RHO
                            if (BDTEST < BDTOL):
                                flag = 99999
                                break

                    if (flag == 99999):
                        flag = 650
                    else:
                        flag = 680
                    break

            else:
                NTRITS = NTRITS + 1
                flag = 90

            #
            #     Severe cancellation is likely to occur if XOPT is too far from XBASE.
            #     If the following test holds, then XBASE is shifted so that XOPT becomes
            #     zero. The appropriate changes are made to BMAT and to the second
            #     derivatives of the current model, beginning with the changes to BMAT
            #     that do not dep on ZMAT. VLAG is used temporarily for working space.
            #
        elif (flag == 90):
            if (DSQ <= 1.0e-3 * XOPTSQ):
                FRACSQ = 0.25e0 * XOPTSQ
                SUMPQ = ZERO
                for K in range(1, NPT+1):
                    SUMPQ = SUMPQ + PQ[K - 1]
                    SUM = -HALF * XOPTSQ
                    for I in range(1, N+1):
                        SUM = SUM + XPT[K - 1][I - 1] * XOPT[I - 1]

                    W[NPT + K - 1] = SUM
                    TEMP = FRACSQ - HALF * SUM
                    for I in range(1, N+1):
                        W[I - 1] = BMAT[K - 1][I - 1]
                        VLAG[I - 1] = SUM * XPT[K - 1][I - 1] + \
                            TEMP * XOPT[I - 1]
                        IP = NPT + I
                        for J in range(1, I+1):
                            BMAT[IP - 1][J - 1] = BMAT[IP - 1][J - 1] + \
                                W[I - 1] * VLAG[J - 1] + VLAG[I - 1] * W[J - 1]

                #
                #     Then the revisions of BMAT that dep on ZMAT are calculated.
                #
                for JJ in range(1, NPTM+1):
                    SUMZ = ZERO
                    SUMW = ZERO
                    for K in range(1, NPT+1):
                        SUMZ = SUMZ + ZMAT[K - 1][JJ - 1]
                        VLAG[K - 1] = W[NPT + K - 1] * ZMAT[K - 1][JJ - 1]
                        SUMW = SUMW + VLAG[K - 1]

                    for J in range(1, N+1):
                        SUM = (FRACSQ * SUMZ - HALF * SUMW) * XOPT[J - 1]
                        for K in range(1, NPT+1):
                            SUM = SUM + VLAG[K - 1] * XPT[K - 1][J - 1]

                        W[J - 1] = SUM
                        for K in range(1, NPT+1):
                            BMAT[K - 1][J - 1] = BMAT[K - 1][J - 1] + \
                                SUM * ZMAT[K - 1][JJ - 1]

                    for I in range(1, N+1):
                        IP = I + NPT
                        TEMP = W[I - 1]
                        for J in range(1, I+1):
                            BMAT[IP - 1][J - 1] = BMAT[IP -
                                                       1][J - 1] + TEMP * W[J - 1]

                #
                #     The following instructions complete the shift, including the changes
                #     to the second derivative parameters of the quadratic model.
                #
                IH = 0
                for J in range(1, N+1):
                    W[J - 1] = -HALF * SUMPQ * XOPT[J - 1]
                    for K in range(1, NPT+1):
                        W[J - 1] = W[J - 1] + PQ[K - 1] * XPT[K - 1][J - 1]
                        XPT[K - 1][J - 1] = XPT[K - 1][J - 1] - XOPT[J - 1]

                    for I in range(1, J+1):
                        IH = IH + 1
                        HQ[IH - 1] = HQ[IH - 1] + W[I - 1] * \
                            XOPT[J - 1] + XOPT[I - 1] * W[J - 1]
                        BMAT[NPT + I - 1][J - 1] = BMAT[NPT + J - 1][I - 1]

                for I in range(1, N+1):
                    XBASE[I - 1] = XBASE[I - 1] + XOPT[I - 1]
                    XNEW[I - 1] = XNEW[I - 1] - XOPT[I - 1]
                    SL[I - 1] = SL[I - 1] - XOPT[I - 1]
                    SU[I - 1] = SU[I - 1] - XOPT[I - 1]
                    XOPT[I - 1] = ZERO

                XOPTSQ = ZERO

            if (NTRITS == 0):
                flag = 210
            else:
                flag = 230

            #
            #     XBASE is also moved to XOPT by a call of RESCUE. This calculation is

            #     more expensive than the previous shift, because new matrices BMAT and
            #     ZMAT are generated from scratch, which may include the replacement of
            #     interpolation points whose positions seem to be causing near linear
            #     depence in the interpolation conditions. Therefore RESCUE is called
            #     only if rounding errors have reduced by at least a factor of two the
            #     denominator of the formula for updating the H matrix. It provides a
            #     useful safeguard, but is not invoked in most applications of BOBYQA.
            #
        elif (flag == 190):
            NFSAV = NF
            KBASE = KOPT
            (XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU, NF, KOPT,
             VLAG, W[:N], W[N:2 * N], W[2*N:NDIM + N], W[NDIM + N:3 * NDIM]) = \
                rescue(N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
                       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, NF, DELTA, KOPT, VLAG,
                       W[:N], W[N:2 * N], W[2*N:NDIM + N], W[NDIM + N:3 * NDIM])
            #
            #     XOPT is updated now in case the branch below to label 720 is taken.
            #     Any updating of GOPT occurs after the branch below to label 20, which
            #     leads to a trust region iteration as does the branch to label 60.
            #
            XOPTSQ = ZERO
            if (KOPT != KBASE):
                for I in range(1, N+1):
                    XOPT[I - 1] = XPT[KOPT - 1][I - 1]
                    XOPTSQ = XOPTSQ + XOPT[I - 1] ** 2

            if (NF < 0):
                NF = MAXFUN
                if (IPRINT > 0):
                    print(
                        'Return from BOBYQA because CALFUN has been called' +
                        ' MAXFUN times.')

                break

            NRESC = NF
            if (NFSAV < NF):
                NFSAV = NF
                flag = 20
            elif (NTRITS > 0):
                flag = 60
            else:
                flag = 210

            #
            #     Pick two alternative vectors of variables, relative to XBASE, that
            #     are suitable as new positions of the KNEW-th interpolation point.
            #     Firstly, XNEW is set to the point on a line through XOPT and another
            #     interpolation point that minimizes the predicted value of the next
            #     denominator, subject to andXNEW - XOPTand .LEQ. ADELT and to the SL
            #     and SU bounds. Secondly, XALT is set to the best feasible point on
            #     a constrained version of the Cauchy step of the KNEW-th Lagrange
            #     function, the corresponding value of the square of this function
            #     being returned in CAUCHY. The choice between these alternatives is
            #     going to be made when the denominator is calculated.
            #
        elif (flag == 210):
            (XNEW, XALT, ALPHA, CAUCHY, W[:  N], W[N:  NDIM], W[NDIM:  NDIM + 2 * N]) \
                = altmov(N, NPT, XPT, XOPT, BMAT, ZMAT, SL, SU, KOPT, KNEW,
                         ADELT, XNEW, XALT, CAUCHY, W[:  N], W[N:  NDIM], W[NDIM:  NDIM + 2 * N])
            for I in range(1, N+1):
                D[I - 1] = XNEW[I - 1] - XOPT[I - 1]

            flag = 230
            #
            #     Calculate VLAG and BETA for the current choice of D. The scalar
            #     product of D with XPT(K,.) is going to be held in W(NPT+K) for
            #     use when VQUAD is calculated.
            #
        elif (flag == 230):
            for K in range(1, NPT+1):
                SUMA = ZERO
                SUMB = ZERO
                SUM = ZERO
                for J in range(1, N+1):
                    SUMA = SUMA + XPT[K - 1][J - 1] * D[J - 1]
                    SUMB = SUMB + XPT[K - 1][J - 1] * XOPT[J - 1]
                    SUM = SUM + BMAT[K - 1][J - 1] * D[J - 1]

                W[K - 1] = SUMA * (HALF * SUMA + SUMB)
                VLAG[K - 1] = SUM
                W[NPT + K - 1] = SUMA

            BETA = ZERO
            for JJ in range(1, NPTM+1):
                SUM = ZERO
                for K in range(1, NPT+1):
                    SUM = SUM + ZMAT[K - 1][JJ - 1] * W[K - 1]

                BETA = BETA - SUM * SUM
                for K in range(1, NPT+1):
                    VLAG[K - 1] = VLAG[K - 1] + SUM * ZMAT[K - 1][JJ - 1]

            DSQ = ZERO
            BSUM = ZERO
            DX = ZERO
            for J in range(1, N+1):
                DSQ = DSQ + D[J - 1] ** 2
                SUM = ZERO
                for K in range(1, NPT+1):
                    SUM = SUM + W[K - 1] * BMAT[K - 1][J - 1]

                BSUM = BSUM + SUM * D[J - 1]
                JP = NPT + J
                for I in range(1, N+1):
                    SUM = SUM + BMAT[JP - 1][I - 1] * D[I - 1]

                VLAG[JP - 1] = SUM
                BSUM = BSUM + SUM * D[J - 1]
                DX = DX + D[J - 1] * XOPT[J - 1]

            BETA = DX * DX + DSQ * \
                (XOPTSQ + DX + DX + HALF * DSQ) + BETA - BSUM
            VLAG[KOPT - 1] = VLAG[KOPT - 1] + ONE
            #
            #     If NTRITS is zero, the denominator may be increased by replacing
            #     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
            #     rounding errors have damaged the chosen denominator.
            #
            if (NTRITS == 0):
                DENOM = VLAG[KNEW - 1] ** 2 + ALPHA * BETA
                if (DENOM < CAUCHY and CAUCHY > ZERO):
                    for I in range(1, N+1):
                        XNEW[I - 1] = XALT[I - 1]
                        D[I - 1] = XNEW[I - 1] - XOPT[I - 1]

                    CAUCHY = ZERO
                    flag = 230
                elif (DENOM <= HALF * VLAG[KNEW - 1] ** 2):
                    if (NF > NRESC):
                        flag = 190
                    else:
                        if (IPRINT > 0):
                            print('Return from BOBYQA because of much' +
                                  ' cancellation in a denominator.')

                        break

                else:
                    flag = 360

                #
                #     Alternatively, if NTRITS is positive, then set KNEW to the index of
                #     the next interpolation point to be deleted to make room for a trust
                #     \region step. Again RESCUE may be called if rounding errors have damaged
                #     the chosen denominator, which is the reason for attempting to select
                #     KNEW before calculating the next value of the objective function.
                #
            else:
                DELSQ = DELTA * DELTA
                SCADEN = ZERO
                BIGLSQ = ZERO
                KNEW = 0
                for K in range(1, NPT+1):
                    if (K == KOPT):
                        continue

                    HDIAG = ZERO
                    for JJ in range(1, NPTM+1):
                        HDIAG = HDIAG + ZMAT[K - 1][JJ - 1] ** 2

                    DEN = BETA * HDIAG + VLAG[K - 1] ** 2
                    DISTSQ = ZERO
                    for J in range(1, N+1):
                        DISTSQ = DISTSQ + \
                            (XPT[K - 1][J - 1] - XOPT[J - 1]) ** 2

                    TEMP = max(ONE, (DISTSQ / DELSQ) ** 2)
                    if (TEMP * DEN > SCADEN):
                        SCADEN = TEMP * DEN
                        KNEW = K
                        DENOM = DEN

                    BIGLSQ = max(BIGLSQ, TEMP * VLAG[K - 1] ** 2)

                if (SCADEN <= HALF * BIGLSQ):
                    if (NF > NRESC):
                        flag = 190
                    else:
                        if (IPRINT > 0):
                            print('Return from BOBYQA because of much' +
                                  ' cancellation in a denominator.')

                        break

                else:
                    flag = 360

            #
            #     Put the variables for the next calculation of the objective function
            #       in XNEW, with any adjustments for the bounds.
            #
            #
            #     Calculate the value of the objective function at XBASE+XNEW, unless
            #       the limit on the number of calculations of F has been reached.
            #
        elif (flag == 360):
            for I in range(1, N+1):
                X[I - 1] = min(max(XL[I - 1], XBASE[I - 1] +
                               XNEW[I - 1]), XU[I - 1])
                if (XNEW[I - 1] == SL[I - 1]):
                    X[I - 1] = XL[I - 1]

                if (XNEW[I - 1] == SU[I - 1]):
                    X[I - 1] = XU[I - 1]

            if (NF >= MAXFUN):
                if (IPRINT > 0):
                    print('Return from BOBYQA because CALFUN has been' +
                          ' called MAXFUN times.')

                break

            NF = NF + 1
            F = calfun(N, X)
            if (IPRINT == 3):
                print(f'Function number {NF}    F ={F}    ' +
                      f'The corresponding X is: {X[: N]}.')

            if (NTRITS == -1):
                FSAVE = F
                break

            #
            #     Use the quadratic model to predict the change in F due to the step D,
            #       and set DIFF to the error of this prediction.
            #
            FOPT = FVAL[KOPT - 1]
            VQUAD = ZERO
            IH = 0
            for J in range(1, N+1):
                VQUAD = VQUAD + D[J - 1] * GOPT[J - 1]
                for I in range(1, J+1):
                    IH = IH + 1
                    TEMP = D[I - 1] * D[J - 1]
                    if (I == J):
                        TEMP = HALF * TEMP

                    VQUAD = VQUAD + HQ[IH - 1] * TEMP

            for K in range(1, NPT+1):
                VQUAD = VQUAD + HALF * PQ[K - 1] * W[NPT + K - 1] ** 2

            DIFF = F - FOPT - VQUAD
            DIFFC = DIFFB
            DIFFB = DIFFA
            DIFFA = abs(DIFF)
            if (DNORM > RHO):
                NFSAV = NF

            #
            #     Pick the next value of DELTA after a trust region step.
            #
            if (NTRITS > 0):
                if (VQUAD >= ZERO):
                    if (IPRINT > 0):
                        print('Return from BOBYQA because a trust' +
                              ' region step has failed to reduce Q.')

                    break

                RATIO = (F - FOPT) / VQUAD
                if (RATIO <= TENTH):
                    DELTA = min(HALF * DELTA, DNORM)
                elif (RATIO <= 0.7e0):
                    DELTA = max(HALF * DELTA, DNORM)
                else:
                    DELTA = max(HALF * DELTA, DNORM + DNORM)

                if (DELTA <= 1.5e0 * RHO):
                    DELTA = RHO

                #
                #     Recalculate KNEW and DENOM if the new F is less than FOPT.
                #
                if (F < FOPT):
                    KSAV = KNEW
                    DENSAV = DENOM
                    DELSQ = DELTA * DELTA
                    SCADEN = ZERO
                    BIGLSQ = ZERO
                    KNEW = 0
                    for K in range(1, NPT+1):
                        HDIAG = ZERO
                        for JJ in range(1, NPTM+1):
                            HDIAG = HDIAG + ZMAT[K - 1][JJ - 1] ** 2

                        DEN = BETA * HDIAG + VLAG[K - 1] ** 2
                        DISTSQ = ZERO
                        for J in range(1, N+1):
                            DISTSQ = DISTSQ + \
                                (XPT[K - 1][J - 1] - XNEW[J - 1]) ** 2

                        TEMP = max(ONE, (DISTSQ / DELSQ) ** 2)
                        if (TEMP * DEN > SCADEN):
                            SCADEN = TEMP * DEN
                            KNEW = K
                            DENOM = DEN

                        BIGLSQ = max(BIGLSQ, TEMP * VLAG[K - 1] ** 2)

                    if (SCADEN <= HALF * BIGLSQ):
                        KNEW = KSAV
                        DENOM = DENSAV

            #
            #     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
            #     moved. Also update the second derivative terms of the model.
            #
            (BMAT, ZMAT, VLAG, W) = update(
                N, NPT, BMAT, ZMAT, VLAG, BETA, DENOM, KNEW, W)
            IH = 0
            PQOLD = PQ[KNEW - 1]
            PQ[KNEW - 1] = ZERO
            for I in range(1, N+1):
                TEMP = PQOLD * XPT[KNEW - 1][I - 1]
                for J in range(1, I+1):
                    IH = IH + 1
                    HQ[IH - 1] = HQ[IH - 1] + TEMP * XPT[KNEW - 1][J - 1]

            for JJ in range(1, NPTM+1):
                TEMP = DIFF * ZMAT[KNEW - 1][JJ - 1]

                for K in range(1, NPT+1):
                    PQ[K - 1] = PQ[K - 1] + TEMP * ZMAT[K - 1][JJ - 1]

            #
            #     Include the new interpolation point, and make the changes to GOPT at
            #     the old XOPT that are caused by the updating of the quadratic model.
            #
            FVAL[KNEW - 1] = F
            for I in range(1, N+1):
                XPT[KNEW - 1][I - 1] = XNEW[I - 1]
                W[I - 1] = BMAT[KNEW - 1][I - 1]

            for K in range(1, NPT+1):
                SUMA = ZERO

                for JJ in range(1, NPTM+1):
                    SUMA = SUMA + ZMAT[KNEW - 1][JJ - 1] * ZMAT[K - 1][JJ - 1]

                SUMB = ZERO
                for J in range(1, N+1):
                    SUMB = SUMB + XPT[K - 1][J - 1] * XOPT[J - 1]

                TEMP = SUMA * SUMB
                for I in range(1, N+1):
                    W[I - 1] = W[I - 1] + TEMP * XPT[K - 1][I - 1]

            for I in range(1, N+1):
                GOPT[I - 1] = GOPT[I - 1] + DIFF * W[I - 1]

            #
            #     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
            #
            if (F < FOPT):
                KOPT = KNEW
                XOPTSQ = ZERO
                IH = 0
                for J in range(1, N+1):
                    XOPT[J - 1] = XNEW[J - 1]
                    XOPTSQ = XOPTSQ + XOPT[J - 1] ** 2
                    for I in range(1, J+1):
                        IH = IH + 1
                        if (I < J):
                            GOPT[J - 1] = GOPT[J - 1] + HQ[IH - 1] * D[I - 1]

                        GOPT[I - 1] = GOPT[I - 1] + HQ[IH - 1] * D[J - 1]

                for K in range(1, NPT+1):
                    TEMP = ZERO
                    for J in range(1, N+1):
                        TEMP = TEMP + XPT[K - 1][J - 1] * D[J - 1]

                    TEMP = PQ[K - 1] * TEMP
                    for I in range(1, N+1):
                        GOPT[I - 1] = GOPT[I - 1] + TEMP * XPT[K - 1][I - 1]

            #
            #     Calculate the parameters of the least Frobenius norm interpolant to
            #     the current data, the gradient of this interpolant at XOPT being put
            #     into VLAG(NPT+I), I=1,2,...,N.
            #
            if (NTRITS > 0):
                for K in range(1, NPT+1):
                    VLAG[K - 1] = FVAL[K - 1] - FVAL[KOPT - 1]
                    W[K - 1] = ZERO

                for J in range(1, NPTM+1):
                    SUM = ZERO
                    for K in range(1, NPT+1):
                        SUM = SUM + ZMAT[K - 1][J - 1] * VLAG[K - 1]

                    for K in range(1, NPT+1):
                        W[K - 1] = W[K - 1] + SUM * ZMAT[K - 1][J - 1]

                for K in range(1, NPT+1):
                    SUM = ZERO
                    for J in range(1, N+1):
                        SUM = SUM + XPT[K - 1][J - 1] * XOPT[J - 1]

                    W[K + NPT - 1] = W[K - 1]
                    W[K - 1] = SUM * W[K - 1]

                GQSQ = ZERO
                GISQ = ZERO
                for I in range(1, N+1):
                    SUM = ZERO
                    for K in range(1, NPT+1):
                        SUM = SUM + BMAT[K - 1][I - 1] * \
                            VLAG[K - 1] + XPT[K - 1][I - 1] * W[K - 1]

                    if (XOPT[I - 1] == SL[I - 1]):
                        GQSQ = GQSQ + min(ZERO, GOPT[I - 1]) ** 2
                        GISQ = GISQ + min(ZERO, SUM) ** 2
                    elif (XOPT[I - 1] == SU[I - 1]):
                        GQSQ = GQSQ + max(ZERO, GOPT[I - 1]) ** 2
                        GISQ = GISQ + max(ZERO, SUM) ** 2
                    else:
                        GQSQ = GQSQ + GOPT[I - 1] ** 2
                        GISQ = GISQ + SUM * SUM

                    VLAG[NPT + I - 1] = SUM

                #
                #     Test whether to replace the new quadratic model by the least Frobenius
                #     norm interpolant, making the replacement if the test is satisfied.
                #
                ITEST = ITEST + 1
                if (GQSQ < TEN * GISQ):
                    ITEST = 0

                if (ITEST >= 3):
                    for I in range(1, trunc(max(NPT, NH))+1):
                        if (I <= N):
                            GOPT[I - 1] = VLAG[NPT + I - 1]

                        if (I <= NPT):
                            PQ[I - 1] = W(NPT + I)

                        if (I <= NH):
                            HQ[I - 1] = ZERO

                        ITEST = 0

            #
            #     If a trust region step has provided a sufficient decrease in F, then
            #     branch for another trust region calculation. The case NTRITS=0 occurs
            #     when the new interpolation point was reached by an alternative step.
            #
            if (NTRITS == 0 or F <= FOPT + TENTH * VQUAD):
                flag = 60
            else:
                #
                #     Alternatively, find out if the interpolation points are close enough
                #       to the best point so far.
                #
                DISTSQ = max((TWO * DELTA) ** 2, (TEN * RHO) ** 2)
                flag = 650

        elif (flag == 650):
            KNEW = 0
            for K in range(1, NPT+1):
                SUM = ZERO
                for J in range(1, N+1):
                    SUM = SUM + (XPT[K - 1][J - 1] - XOPT[J - 1]) ** 2

                if (SUM > DISTSQ):
                    KNEW = K
                    DISTSQ = SUM

            #
            #     If KNEW is positive, then ALTMOV finds alternative new positions for
            #     the KNEW-th interpolation point within distance ADELT of XOPT. It is
            #     reached via label 90. Otherwise, there is a branch to label 60 for
            #     another trust region iteration, unless the calculations with the
            #     current RHO are complete.
            #
            if (KNEW > 0):
                DIST = sqrt(DISTSQ)
                if (NTRITS == -1):
                    DELTA = min(TENTH * DELTA, HALF * DIST)
                    if (DELTA <= 1.5e0 * RHO):
                        DELTA = RHO

                NTRITS = 0
                ADELT = max(min(TENTH * DIST, DELTA), RHO)
                DSQ = ADELT * ADELT
                flag = 90
            elif (NTRITS == -1):
                flag = 680
            elif (RATIO > ZERO or max(DELTA, DNORM) > RHO):
                flag = 60
            else:
                flag = 680

            #
            #     The calculations with the current value of RHO are complete. Pick the
            #       next values of RHO and DELTA.
            #
        elif (flag == 680):
            if (RHO > RHOEND):
                DELTA = HALF * RHO
                RATIO = RHO / RHOEND
                if (RATIO <= 16.0e0):
                    RHO = RHOEND
                elif (RATIO <= 250.0e0):
                    RHO = sqrt(RATIO) * RHOEND
                else:
                    RHO = TENTH * RHO

                DELTA = max(DELTA, RHO)
                if (IPRINT >= 2):
                    print(f'New RHO = {RHO} Number of function values = {NF}.')
                    print(f'Least value of F = {FVAL[KOPT - 1]}.')
                    print(
                        f'The corresponding X is: {[XBASE[I] + XOPT[I] for I in range(N)]}')

                NTRITS = 0
                NFSAV = NF
                flag = 60
                #
                #     Return from the calculation, after another Newton-Raphson step, if
                #       it is too short to have been tried before.
                #
            elif (NTRITS == -1):
                flag = 360
            else:
                break

    if (FVAL[KOPT - 1] <= FSAVE):
        for I in range(1, N+1):
            X[I - 1] = min(max(XL[I - 1], XBASE[I - 1] +
                           XOPT[I - 1]), XU[I - 1])
            if (XOPT[I - 1] == SL[I - 1]):
                X[I - 1] = XL[I - 1]

            if (XOPT[I - 1] == SU[I - 1]):
                X[I - 1] = XU[I - 1]
        F = FVAL[KOPT - 1]

    if (IPRINT >= 1):
        print(f'At the return from BOBYQA. Number of function values = {NF}.')
        print(f'Least value of F = {F}.')
        print(f'The corresponding X is: {X[: N]}.')


def update(N, NPT, BMAT, ZMAT, VLAG, BETA, DENOM, KNEW, W):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
    #
    #     The arrays BMAT and ZMAT are updated, as required by the new position
    #     of the interpolation point that has the index KNEW. The vector VLAG has
    #     N+NPT components, set on entry to the first NPT and last N components
    #     of the product Hw in equation (4.11) of the Powell (2006) paper on
    #     NEWUOA. Further, BETA is set on entry to the value of the parameter
    #     with that name, and DENOM is set to the denominator of the updating
    #     formula. Elements of ZMAT may be treated as zero if their moduli are
    #     at most ZTEST. The first NDIM elements of W are used for working space.
    #
    #     Set some constants.
    #
    ONE = 1.0e0
    ZERO = 0.0e0
    NPTM = NPT - N - 1
    ZTEST = ZERO
    for K in range(1, NPT):
        for J in range(1, NPTM+1):
            ZTEST = max(ZTEST, abs(ZMAT[K - 1][J - 1]))

    ZTEST = 1.0e-20 * ZTEST
    #
    #     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
    #
    for J in range(2, NPTM+1):
        if (abs(ZMAT[KNEW - 1][J - 1]) > ZTEST):
            TEMP = sqrt(ZMAT[KNEW - 1][0] ** 2 + ZMAT[KNEW - 1][J - 1] ** 2)
            TEMPA = ZMAT[KNEW - 1][0] / TEMP
            TEMPB = ZMAT[KNEW - 1][J - 1] / TEMP
            for I in range(1, NPT+1):
                TEMP = TEMPA * ZMAT[I - 1][0] + TEMPB * ZMAT[I - 1][J - 1]
                ZMAT[I - 1][J - 1] = TEMPA * ZMAT[I - 1][J - 1] - \
                    TEMPB * ZMAT[I - 1][0]
                ZMAT[I - 1][0] = TEMP

        ZMAT[KNEW - 1][J - 1] = ZERO

    #
    #     Put the first NPT components of the KNEW-th column of HLAG into W,
    #     and calculate the parameters of the updating formula.
    #
    for I in range(1, NPT+1):
        W[I - 1] = ZMAT[KNEW - 1][0] * ZMAT[I - 1][0]

    ALPHA = W[KNEW - 1]
    TAU = VLAG[KNEW - 1]
    VLAG[KNEW - 1] = VLAG[KNEW - 1] - ONE
    #
    #     Complete the updating of ZMAT.
    #
    TEMP = sqrt(DENOM)
    TEMPB = ZMAT[KNEW - 1][0] / TEMP
    TEMPA = TAU / TEMP
    for I in range(1, NPT+1):
        ZMAT[I - 1][0] = TEMPA * ZMAT[I - 1][0] - TEMPB * VLAG[I - 1]

    #
    #     Finally, update the matrix BMAT.
    #
    for J in range(1, N+1):
        JP = NPT + J
        W[JP - 1] = BMAT[KNEW - 1][J - 1]
        TEMPA = (ALPHA * VLAG[JP - 1] - TAU * W[JP - 1]) / DENOM
        TEMPB = (-BETA * W[JP - 1] - TAU * VLAG[JP - 1]) / DENOM
        for I in range(1, JP+1):
            BMAT[I - 1][J - 1] = BMAT[I - 1][J - 1] + \
                TEMPA * VLAG[I - 1] + TEMPB * W[I - 1]
            if (I > NPT):
                BMAT[JP - 1][I - NPT - 1] = BMAT[I - 1][J - 1]

    return (BMAT, ZMAT, VLAG, W)


def trsbox(N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL, SU, DELTA,
           XNEW, D, GNEW, XBDI, S, HS, HRED):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*),
    #      1     XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
    #
    #     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
    #       meanings as the corresponding arguments of BOBYQB.
    #     DELTA is the trust region radius for the present calculation, which
    #       seeks a small value of the quadratic model within distance DELTA of
    #       XOPT subject to the bounds on the variables.
    #     XNEW will be set to a new vector of variables that is approximately
    #       the one that minimizes the quadratic model within the trust region
    #       subject to the SL and SU constraints on the variables. It satisfies
    #       as equations the bounds that become active during the calculation.
    #     D is the calculated trial step from XOPT, generated iteratively from an
    #       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
    #     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
    #       when D is updated.
    #     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
    #       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
    #       I-th variable has become fixed at a bound, the bound being SL(I) or
    #       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
    #       information is accumulated during the construction of XNEW.
    #     The arrays S, HS and HRED are also used for working space. They hold the
    #       current search direction, and the changes in the gradient of Q along S
    #       and the reduced D, respectively, where the reduced D is the same as D,
    #       except that the components of the fixed variables are zero.
    #     DSQ will be set to the square of the length of XNEW-XOPT.
    #     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
    #       it is set to the least curvature of H that occurs in the conjugate
    #       gradient searches that are not restricted by any constraints. The
    #       value CRVMIN=-1.0e0 is set, however, if all of these searches are
    #       constrained.
    #
    #     A version of the truncated conjugate gradient is applied. If a line
    #     search is restricted by a constraint, then the procedure is restarted,
    #     the values of the variables that are at their bounds being fixed. If
    #     the trust region boundary is reached, then further changes may be made
    #     to D, each one being in the two dimensional space that is spanned
    #     by the current D and the gradient of Q at XOPT+D, staying on the trust
    #     \region boundary. Termination occurs when the reduction in Q seems to
    #     be close to the greatest reduction that can be achieved.
    #
    #     Set some constants.
    #
    HALF = 0.5e0
    ONE = 1.0e0
    ONEMIN = -1.0e0
    ZERO = 0.0e0
    #
    #     The sign of GOPT(I) gives the sign of the change to the I-th variable
    #     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
    #     or not to fix the I-th variable at one of its bounds initially, with
    #     NACT being set to the number of fixed variables. D and GNEW are also
    #     set for the first iteration. DELSQ is the upper bound on the sum of
    #     squares of the free variables. QRED is the reduction in Q so far.
    #
    ITERC = 0
    NACT = 0
    for I in range(1, N+1):
        XBDI[I - 1] = ZERO
        if (XOPT[I-1] <= SL[I-1]):
            if (GOPT[I - 1] >= ZERO):
                XBDI[I - 1] = ONEMIN

        elif (XOPT[I - 1] >= SU[I - 1]):
            if (GOPT[I - 1] <= ZERO):
                XBDI[I - 1] = ONE

        if (XBDI[I-1] != ZERO):
            NACT = NACT + 1

        D[I - 1] = ZERO
        GNEW[I - 1] = GOPT[I - 1]

    DELSQ = DELTA * DELTA
    QRED = ZERO
    CRVMIN = ONEMIN
    #
    #     Set the next search direction of the conjugate gradient method. It is
    #     the steepest descent direction initially and when the iterations are
    #     restarted because a variable has just been fixed by a bound, and of
    #     course the components of the fixed variables are zero. ITERMAX is an
    #     upper bound on the indices of the conjugate gradient iterations.
    #
    flag = 20
    while (True):
        if (flag == 20):
            BETA = ZERO
            flag = 30
        elif (flag == 30):
            STEPSQ = ZERO
            for I in range(1, N+1):
                if (XBDI[I - 1] != ZERO):
                    S[I - 1] = ZERO
                elif (BETA == ZERO):
                    S[I - 1] = -GNEW[I - 1]
                else:
                    S[I - 1] = BETA * S[I - 1] - GNEW[I - 1]

                STEPSQ = STEPSQ + S[I - 1] ** 2

            if (STEPSQ == ZERO):
                break

            if (BETA == ZERO):
                GREDSQ = STEPSQ
                ITERMAX = ITERC + N - NACT

            if (GREDSQ * DELSQ <= 1.0e-4 * QRED * QRED):
                break

            #
            #     Multiply the search direction by the second derivative matrix of Q and
            #     calculate some scalars for the choice of steplength. Then set BLEN to
            #     the length of the the step to the trust region boundary and STPLEN to
            #     the steplength, ignoring the simple bounds.
            #
            flag = 210
        elif (flag == 50):
            RESID = DELSQ
            DS = ZERO
            SHS = ZERO
            for I in range(1, N+1):
                if (XBDI[I - 1] == ZERO):
                    RESID = RESID - D[I - 1] ** 2
                    DS = DS + S[I - 1] * D[I - 1]
                    SHS = SHS + S[I - 1] * HS[I - 1]

            if (RESID <= ZERO):
                flag = 90
            else:
                TEMP = sqrt(STEPSQ * RESID + DS * DS)
                if (DS < ZERO):
                    BLEN = (TEMP - DS) / STEPSQ
                else:
                    BLEN = RESID / (TEMP + DS)

                STPLEN = BLEN
                if (SHS > ZERO):
                    STPLEN = min(BLEN, GREDSQ / SHS)

                #
                #     Reduce STPLEN if necessary in order to preserve the simple bounds,
                #     letting IACT be the index of the new constrained variable.
                #
                IACT = 0
                for I in range(1, N+1):
                    if (S[I - 1] != ZERO):
                        XSUM = XOPT[I - 1] + D[I - 1]
                        if (S[I - 1] > ZERO):
                            TEMP = (SU[I - 1] - XSUM) / S[I - 1]
                        else:
                            TEMP = (SL[I - 1] - XSUM) / S[I - 1]

                        if (TEMP < STPLEN):
                            STPLEN = TEMP
                            IACT = I

                #
                #     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
                #
                SDEC = ZERO
                if (STPLEN > ZERO):
                    ITERC = ITERC + 1
                    TEMP = SHS / STEPSQ
                    if (IACT == 0 and TEMP > ZERO):
                        CRVMIN = min(CRVMIN, TEMP)
                        if (CRVMIN == ONEMIN):
                            CRVMIN = TEMP

                    GGSAV = GREDSQ
                    GREDSQ = ZERO
                    for I in range(1, N+1):
                        GNEW[I - 1] = GNEW[I - 1] + STPLEN * HS[I - 1]
                        if (XBDI[I - 1] == ZERO):
                            GREDSQ = GREDSQ + GNEW[I - 1] ** 2

                        D[I - 1] = D[I - 1] + STPLEN * S[I - 1]

                    SDEC = max(STPLEN * (GGSAV - HALF * STPLEN * SHS), ZERO)
                    QRED = QRED + SDEC

                #
                #     Restart the conjugate gradient method if it has hit a new bound.
                #
                if (IACT > 0):
                    NACT = NACT + 1
                    XBDI[IACT - 1] = ONE
                    if (S[IACT - 1] < ZERO):
                        XBDI[IACT - 1] = ONEMIN

                    DELSQ = DELSQ - D[IACT - 1] ** 2
                    if (DELSQ <= ZERO):
                        flag = 90
                    else:
                        flag = 20

                else:
                    #
                    #     If STPLEN is less than BLEN, then either apply another conjugate
                    #     gradient iteration or RETURN.
                    #
                    if (STPLEN < BLEN):
                        if (ITERC == ITERMAX or SDEC <= 0.01e0 * QRED):
                            break

                        BETA = GREDSQ / GGSAV
                        flag = 30
                    else:
                        flag = 90

        elif (flag == 90):
            CRVMIN = ZERO
            flag = 100
            #
            #     Prepare for the alternative iteration by calculating some scalars and
            #     by multiplying the reduced D by the second derivative matrix of Q.
            #
        elif (flag == 100):
            if (NACT >= N - 1):
                break

            DREDSQ = ZERO
            DREDG = ZERO
            GREDSQ = ZERO
            for I in range(1, N+1):
                if (XBDI[I - 1] == ZERO):
                    DREDSQ = DREDSQ + D[I - 1] ** 2
                    DREDG = DREDG + D[I - 1] * GNEW[I - 1]
                    GREDSQ = GREDSQ + GNEW[I - 1] ** 2
                    S[I - 1] = D[I - 1]
                else:
                    S[I - 1] = ZERO

            ITCSAV = ITERC
            flag = 210
            #
            #     The following instructions multiply the current S-vector by the second
            #     derivative matrix of the quadratic model, putting the product in HS.
            #     They are reached from three different parts of the software above and
            #     they can be regarded as an external subroutine.
            #
        elif (flag == 210):
            IH = 0
            for J in range(1, N+1):
                HS[J - 1] = ZERO
                for I in range(1, J+1):
                    IH = IH + 1
                    if (I < J):
                        HS[J - 1] = HS[J - 1] + HQ[IH - 1] * S[I - 1]

                    HS[I - 1] = HS[I - 1] + HQ[IH - 1] * S[J - 1]

            for K in range(1, NPT+1):
                if (PQ[K - 1] != ZERO):
                    TEMP = ZERO
                    for J in range(1, N+1):
                        TEMP = TEMP + XPT[K - 1][J - 1] * S[J - 1]

                    TEMP = TEMP * PQ[K - 1]
                    for I in range(1, N+1):
                        HS[I - 1] = HS[I - 1] + TEMP * XPT[K - 1][I - 1]

            if (CRVMIN != ZERO):
                flag = 50
            elif (ITERC > ITCSAV):
                flag = 150
            else:
                for I in range(1, N+1):
                    HRED[I - 1] = HS[I - 1]

                flag = 120

            #
            #     Let the search direction S be a linear combination of the reduced D
            #     and the reduced G that is orthogonal to the reduced D.
            #
        elif (flag == 120):
            ITERC = ITERC + 1
            TEMP = GREDSQ * DREDSQ - DREDG * DREDG
            if (TEMP <= 1.0e-4 * QRED * QRED):
                break

            TEMP = sqrt(TEMP)
            for I in range(1, N+1):
                if (XBDI[I - 1] == ZERO):
                    S[I - 1] = (DREDG * D[I - 1] - DREDSQ * GNEW[I - 1]) / TEMP
                else:
                    S[I - 1] = ZERO

            SREDG = -TEMP
            #
            #     By considering the simple bounds on the variables, calculate an upper
            #     bound on the tangent of half the angle of the alternative iteration,
            #     namely ANGBD, except that, if already a free variable has reached a
            #     bound, there is a branch back to label 100 after fixing that variable.
            #
            ANGBD = ONE
            IACT = 0
            for I in range(1, N+1):
                if (XBDI[I - 1] == ZERO):
                    TEMPA = XOPT[I - 1] + D[I - 1] - SL[I - 1]
                    TEMPB = SU[I - 1] - XOPT[I - 1] - D[I - 1]
                    if (TEMPA <= ZERO):
                        NACT = NACT + 1
                        XBDI[I - 1] = ONEMIN
                        flag = 100
                        break
                    elif (TEMPB <= ZERO):
                        NACT = NACT + 1
                        XBDI[I - 1] = ONE
                        flag = 100
                        break

                    SSQ = D[I - 1] ** 2 + S[I - 1] ** 2
                    TEMP = SSQ - (XOPT[I - 1] - SL[I - 1]) ** 2
                    if (TEMP > ZERO):
                        TEMP = sqrt(TEMP) - S[I - 1]
                        if (ANGBD * TEMP > TEMPA):
                            ANGBD = TEMPA / TEMP
                            IACT = I
                            XSAV = ONEMIN

                    TEMP = SSQ - (SU[I - 1] - XOPT[I - 1]) ** 2
                    if (TEMP > ZERO):
                        TEMP = sqrt(TEMP) + S[I - 1]
                        if (ANGBD * TEMP > TEMPB):
                            ANGBD = TEMPB / TEMP
                            IACT = I
                            XSAV = ONE

                flag = 210

            #
            #     Calculate HHD and some curvatures for the alternative iteration.
            #
        elif (flag == 150):
            SHS = ZERO
            DHS = ZERO
            DHD = ZERO
            for I in range(1, N+1):
                if (XBDI[I - 1] == ZERO):
                    SHS = SHS + S[I - 1] * HS[I - 1]
                    DHS = DHS + D[I - 1] * HS[I - 1]
                    DHD = DHD + D[I - 1] * HRED[I - 1]

            #
            #     Seek the greatest reduction in Q for a range of equally spaced values
            #     of ANGT in [0,ANGBD - 1], where ANGT is the tangent of half the angle of
            #     the alternative iteration.
            #
            REDMAX = ZERO
            ISAV = 0
            REDSAV = ZERO
            IU = trunc(17.0e0 * ANGBD + 3.1e0)
            for I in range(1, IU+1):
                ANGT = ANGBD * I / IU
                STH = (ANGT + ANGT) / (ONE + ANGT * ANGT)
                TEMP = SHS + ANGT * (ANGT * DHD - DHS - DHS)
                REDNEW = STH * (ANGT * DREDG - SREDG - HALF * STH * TEMP)
                if (REDNEW > REDMAX):
                    REDMAX = REDNEW
                    ISAV = I
                    RDPREV = REDSAV
                elif (I == ISAV + 1):
                    RDNEXT = REDNEW

                REDSAV = REDNEW

            #
            #     Return if the reduction is zero. Otherwise, set the sine and cosine
            #     of the angle of the alternative iteration, and calculate SDEC.
            #
            if (ISAV == 0):
                break

            if (ISAV < IU):
                TEMP = (RDNEXT - RDPREV) / (REDMAX + REDMAX - RDPREV - RDNEXT)
                ANGT = ANGBD * (ISAV + HALF * TEMP) / IU

            CTH = (ONE - ANGT * ANGT) / (ONE + ANGT * ANGT)
            STH = (ANGT + ANGT) / (ONE + ANGT * ANGT)
            TEMP = SHS + ANGT * (ANGT * DHD - DHS - DHS)
            SDEC = STH * (ANGT * DREDG - SREDG - HALF * STH * TEMP)
            if (SDEC <= ZERO):
                break

            #
            #     Update GNEW, D and HRED. If the angle of the alternative iteration
            #     is restricted by a bound on a free variable, that variable is fixed
            #     at the bound.
            #
            DREDG = ZERO
            GREDSQ = ZERO
            for I in range(1, N+1):
                GNEW[I - 1] = GNEW[I - 1] + \
                    (CTH - ONE) * HRED[I - 1] + STH * HS[I - 1]
                if (XBDI[I - 1] == ZERO):
                    D[I - 1] = CTH * D[I - 1] + STH * S[I - 1]
                    DREDG = DREDG + D[I - 1] * GNEW[I - 1]
                    GREDSQ = GREDSQ + GNEW[I - 1] ** 2

                HRED[I - 1] = CTH * HRED[I - 1] + STH * HS[I - 1]

            QRED = QRED + SDEC
            if (IACT > 0 and ISAV == IU):
                NACT = NACT + 1
                XBDI[IACT - 1] = XSAV
                flag = 100
                #
                #     If SDEC is sufficiently small, then RETURN after setting XNEW to
                #     XOPT+D, giving careful attention to the bounds.
                #
            elif (SDEC > 0.01e0 * QRED):
                flag = 120
            else:
                break

    DSQ = ZERO
    for I in range(1, N+1):
        XNEW[I - 1] = max(min(XOPT[I - 1] + D[I - 1], SU[I - 1]), SL[I - 1])
        if (XBDI[I - 1] == ONEMIN):
            XNEW[I - 1] = SL[I - 1]

        if (XBDI[I - 1] == ONE):
            XNEW[I - 1] = SU[I - 1]

        D[I - 1] = XNEW[I - 1] - XOPT[I - 1]
        DSQ = DSQ + D[I - 1] ** 2

    return (XNEW, D, GNEW, XBDI, S, HS, HRED, DSQ, CRVMIN)


def rescue(N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT,
           FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL, SU, NF, DELTA,
           KOPT, VLAG, PTSAUX0, PTSAUX1, PTSID, W):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*),
    #      1     GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*),
    #      2     VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
    #
    #     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
    #       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
    #       the corresponding arguments of BOBYQB on the entry to RESCUE.
    #     NF is maintained as the number of calls of CALFUN so far, except that
    #       NF is set to -1 if the value of MAXFUN prevents further progress.
    #     KOPT is maintained so that FVAL(KOPT) is the least calculated function
    #       value. Its correct value must be given on entry. It is updated if a
    #       new least function value is found, but the corresponding changes to
    #       XOPT and GOPT have to be made later by the calling program.
    #     DELTA is the current trust region radius.
    #     VLAG is a working space vector that will be used for the values of the
    #       provisional Lagrange functions at each of the interpolation points.
    #       They are part of a product that requires VLAG to be of length NDIM.
    #     PTSAUX is also a working space array. For J=1:2,...,N, PTSAUX(1,J) and
    #       PTSAUX(2,J) specify the two positions of provisional interpolation
    #       points when a nonzero step is taken along e_J (the J-th coordinate
    #       direction) through XBASE+XOPT, as specified below. Usually these
    #       steps have length DELTA, but other lengths are chosen if necessary
    #       in order to satisfy the given bounds on the variables.
    #     PTSID is also a working space array. It has NPT components that denote
    #       provisional new positions of the original interpolation points, in
    #       case changes are needed to restore the linear independence of the
    #       interpolation conditions. The K-th point is a candidate for change
    #       if and only if PTSID(K) is nonzero. In this case let p and q be the
    #       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
    #       and q are both positive, the step from XBASE+XOPT to the new K-th
    #       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
    #       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
    #       p=0, respectively.
    #     The first NDIM+NPT elements of the array W are used for working space.
    #     The final elements of BMAT and ZMAT are set in a well-conditioned way
    #       to the values that are appropriate for the new interpolation points.
    #     The elements of GOPT, HQ and PQ are also revised to the values that are
    #       appropriate to the final quadratic model.
    #
    #     Set some constants.
    #
    HALF = 0.5e0
    ONE = 1.0e0
    ZERO = 0.0e0
    NP = N + 1
    SFRAC = HALF / (NP)
    NPTM = NPT - NP
    #
    #     Shift the interpolation points so that XOPT becomes the origin, and set
    #     the elements of ZMAT to zero. The value of SUMPQ is required in the
    #     updating of HQ below. The squares of the distances from XOPT to the
    #     other interpolation points are set at the  of W. Increments of WINC
    #     may be added later to these squares to balance the consideration of
    #     the choice of point that is going to become current.
    #
    SUMPQ = ZERO
    WINC = ZERO
    for K in range(1, NPT+1):
        DISTSQ = ZERO
        for J in range(1, N+1):
            XPT[K - 1][J - 1] = XPT[K - 1][J - 1] - XOPT[J - 1]
            DISTSQ = DISTSQ + XPT[K - 1][J - 1] ** 2

        SUMPQ = SUMPQ + PQ[K - 1]
        W[NDIM + K - 1] = DISTSQ
        WINC = max(WINC, DISTSQ)
        for J in range(1, NPTM+1):
            ZMAT[K - 1][J - 1] = ZERO

    #
    #     Update HQ so that HQ and PQ define the second derivatives of the model
    #     after XBASE has been shifted to the trust region centre.
    #
    IH = 0
    for J in range(1, N+1):
        W[J - 1] = HALF * SUMPQ * XOPT[J - 1]
        for K in range(1, NPT+1):
            W[J - 1] = W[J - 1] + PQ[K - 1] * XPT[K - 1][J - 1]

        for I in range(1, J+1):
            IH = IH + 1
            HQ[IH - 1] = HQ[IH - 1] + W[I - 1] * \
                XOPT[J - 1] + W[J - 1] * XOPT[I - 1]

    #
    #     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
    #     also set the elements of PTSAUX.
    #
    for J in range(1, N+1):
        XBASE[J - 1] = XBASE[J - 1] + XOPT[J - 1]
        SL[J - 1] = SL[J - 1] - XOPT[J - 1]
        SU[J - 1] = SU[J - 1] - XOPT[J - 1]
        XOPT[J - 1] = ZERO
        PTSAUX0[J - 1] = min(DELTA, SU[J - 1])
        PTSAUX1[J - 1] = max(-DELTA, SL[J - 1])
        if (PTSAUX0[J - 1] + PTSAUX1[J - 1] < ZERO):
            TEMP = PTSAUX0[J - 1]
            PTSAUX0[J - 1] = PTSAUX1[J - 1]
            PTSAUX1[J - 1] = TEMP

        if (abs(PTSAUX1[J - 1]) < HALF * abs(PTSAUX0[J - 1])):
            PTSAUX1[J - 1] = HALF * PTSAUX0[J - 1]

        for I in range(1, NDIM+1):
            BMAT[I - 1][J - 1] = ZERO

    FBASE = FVAL[KOPT - 1]
    #
    #     Set the identifiers of the artificial interpolation points that are
    #     along a coordinate direction from XOPT, and set the corresponding
    #     nonzero elements of BMAT and ZMAT.
    #
    PTSID[0] = SFRAC
    for J in range(1, N+1):
        JP = J + 1
        JPN = JP + N
        PTSID[JP - 1] = J + SFRAC
        if (JPN <= NPT):
            PTSID[JPN - 1] = J / NP + SFRAC
            TEMP = ONE / (PTSAUX0[J - 1] - PTSAUX1[J - 1])
            BMAT[JP - 1][J - 1] = -TEMP + ONE / PTSAUX0[J - 1]
            BMAT[JPN - 1][J - 1] = TEMP + ONE / PTSAUX1[J - 1]
            BMAT[0][J - 1] = -BMAT[JP - 1][J - 1] - BMAT[JPN - 1][J - 1]
            ZMAT[0][J - 1] = sqrt(2.0e0) / \
                abs(PTSAUX0[J - 1] * PTSAUX1[J - 1])
            ZMAT[JP - 1][J - 1] = ZMAT[0][J - 1] * \
                PTSAUX1[J - 1] * TEMP
            ZMAT[JPN - 1][J - 1] = -ZMAT[0][J - 1] * \
                PTSAUX0[J - 1] * TEMP
        else:
            BMAT[0][J - 1] = -ONE / PTSAUX0[J - 1]
            BMAT[JP - 1][J - 1] = ONE / PTSAUX0[J - 1]
            BMAT[J + NPT - 1][J - 1] = -HALF * PTSAUX0[J - 1] ** 2

    #
    #     Set any remaining identifiers with their nonzero elements of ZMAT.
    #
    if (NPT >= N + NP):
        for K in range((2 * NP), NPT+1):
            IW = trunc(((K - NP) - HALF) / N)
            IP = K - NP - IW * N
            IQ = IP + IW
            if (IQ > N):
                IQ = IQ - N

            PTSID[K - 1] = IP + IQ / NP + SFRAC
            TEMP = ONE / (PTSAUX0[IP - 1] * PTSAUX0[IQ - 1])
            ZMAT[0][K - NP - 1] = TEMP
            ZMAT[IP][K - NP - 1] = -TEMP
            ZMAT[IQ][K - NP - 1] = -TEMP
            ZMAT[K - 1][K - NP - 1] = TEMP

    NREM = NPT
    KOLD = 1
    KNEW = KOPT
    #
    #     Reorder the provisional points in the way that exchanges PTSID(KOLD)
    #     with PTSID(KNEW).
    #
    while (True):
        for J in range(1, N+1):
            TEMP = BMAT[KOLD - 1][J - 1]
            BMAT[KOLD - 1][J - 1] = BMAT[KNEW - 1][J - 1]
            BMAT[KNEW - 1][J - 1] = TEMP

        for J in range(1, NPTM+1):
            TEMP = ZMAT[KOLD - 1][J - 1]
            ZMAT[KOLD - 1][J - 1] = ZMAT[KNEW - 1][J - 1]
            ZMAT[KNEW - 1][J - 1] = TEMP

        PTSID[KOLD - 1] = PTSID[KNEW - 1]
        PTSID[KNEW - 1] = ZERO
        W[NDIM + KNEW - 1] = ZERO
        NREM = NREM - 1
        if (KNEW != KOPT):
            TEMP = VLAG[KOLD - 1]
            VLAG[KOLD - 1] = VLAG[KNEW - 1]
            VLAG[KNEW - 1] = TEMP
            #
            #     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
            #     interpolation point can be changed from provisional to original. The
            #     branch to label 350 occurs if all the original points are reinstated.
            #     The nonnegative values of W(NDIM+K) are required in the search below.
            #
            (BMAT, ZMAT, VLAG, W) = update(
                N, NPT, BMAT, ZMAT, VLAG, BETA, DENOM, KNEW, W)
            if (NREM == 0):
                return (XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL,
                        SU, NF, KOPT, VLAG, PTSAUX0, PTSAUX1, PTSID, W)

            for K in range(1, NPT+1):
                W[NDIM + K - 1] = abs(W[NDIM + K - 1])

        #
        #     Pick the index KNEW of an original interpolation point that has not
        #     yet replaced one of the provisional interpolation points, giving
        #     attention to the closeness to XOPT and to previous tries with KNEW.
        #
        flag = False
        while (True):
            DSQMIN = ZERO
            for K in range(1, NPT+1):
                if (W[NDIM + K - 1] > ZERO):
                    if (DSQMIN == ZERO or W[NDIM + K - 1] < DSQMIN):
                        KNEW = K
                        DSQMIN = W[NDIM + K - 1]

            if (DSQMIN == ZERO):
                flag = True
                break

            #
            #     Form the W-vector of the chosen original interpolation point.
            #
            for J in range(1, N+1):
                W[NPT + J - 1] = XPT[KNEW - 1][J - 1]

            for K in range(1, NPT+1):
                SUM = ZERO
                if (K == KOPT):
                    continue
                elif (PTSID[K - 1] == ZERO):
                    for J in range(1, N+1):
                        SUM = SUM + W[NPT + J - 1] * XPT[K - 1][J - 1]
                else:
                    IP = trunc(PTSID[K - 1])
                    if (IP > 0):
                        SUM = W[NPT + IP - 1] * PTSAUX1[IP - 1]

                    IQ = trunc(NP * PTSID[K - 1] - (IP * NP))
                    if (IQ > 0):
                        IW = 1
                        if (IP == 0):
                            IW = 2
                        if (IW == 1):
                            SUM = SUM + W[NPT + IQ - 1] * PTSAUX0[IQ - 1]
                        elif (IW == 2):
                            SUM = SUM + W[NPT + IQ - 1] * PTSAUX1[IQ - 1]

                W[K - 1] = HALF * SUM * SUM

            #
            #     Calculate VLAG and BETA for the required updating of the H matrix if
            #     XPT(KNEW,.) is reinstated in the set of interpolation points.
            #
            for K in range(1, NPT+1):
                SUM = ZERO
                for J in range(1, N+1):
                    SUM = SUM + BMAT[K - 1][J - 1] * W[NPT + J - 1]

                VLAG[K - 1] = SUM

            BETA = ZERO
            for J in range(1, NPTM+1):
                SUM = ZERO
                for K in range(1, NPT+1):
                    SUM = SUM + ZMAT[K - 1][J - 1] * W[K - 1]

                BETA = BETA - SUM * SUM
                for K in range(1, NPT+1):
                    VLAG[K - 1] = VLAG[K - 1] + SUM * ZMAT[K - 1][J - 1]

            BSUM = ZERO
            DISTSQ = ZERO
            for J in range(1, N+1):
                SUM = ZERO
                for K in range(1, NPT+1):
                    SUM = SUM + BMAT[K - 1][J - 1] * W[K - 1]

                JP = J + NPT
                BSUM = BSUM + SUM * W[JP - 1]
                for IP in range((NPT + 1), NDIM+1):
                    SUM = SUM + BMAT[IP - 1][J - 1] * W[IP - 1]

                BSUM = BSUM + SUM * W[JP - 1]
                VLAG[JP - 1] = SUM
                DISTSQ = DISTSQ + XPT[KNEW - 1][J - 1] ** 2

            BETA = HALF * DISTSQ * DISTSQ + BETA - BSUM
            VLAG[KOPT - 1] = VLAG[KOPT - 1] + ONE
            #
            #     KOLD is set to the index of the provisional interpolation point that is
            #     going to be deleted to make way for the KNEW-th original interpolation
            #     point. The choice of KOLD is governed by the avoidance of a small value
            #     of the denominator in the updating calculation of UPDATE.
            #
            DENOM = ZERO
            VLMXSQ = ZERO
            for K in range(1, NPT+1):
                if (PTSID[K - 1] != ZERO):
                    HDIAG = ZERO
                    for J in range(1, NPTM+1):
                        HDIAG = HDIAG + ZMAT[K - 1][J - 1] ** 2

                    DEN = BETA * HDIAG + VLAG[K - 1] ** 2
                    if (DEN > DENOM):
                        KOLD = K
                        DENOM = DEN

                VLMXSQ = max(VLMXSQ, VLAG[K - 1] ** 2)

            if (DENOM <= 1.0e-2 * VLMXSQ):
                W[NDIM + KNEW - 1] = -W[NDIM + KNEW - 1] - WINC
            else:
                break

        if flag:
            break

    #
    #     When label 260 is reached, all the final positions of the interpolation
    #     points have been chosen although any changes have not been included yet
    #     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
    #     from the shift of XBASE, the updating of the quadratic model remains to
    #     be done. The following cycle through the new interpolation points begins
    #     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
    #     except that a return occurs if MAXFUN prohibits another value of F.
    #
    for KPT in range(1, NPT+1):
        if (PTSID[KPT - 1] == ZERO):
            continue

        if (NF >= MAXFUN):
            NF = -1
            return (XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL,
                    SU, NF, KOPT, VLAG, PTSAUX0, PTSAUX1, PTSID, W)

        IH = 0
        for J in range(1, N+1):
            W[J - 1] = XPT[KPT - 1][J - 1]
            XPT[KPT - 1][J - 1] = ZERO
            TEMP = PQ[KPT - 1] * W[J - 1]
            for I in range(1, J+1):
                IH = IH + 1
                HQ[IH - 1] = HQ[IH - 1] + TEMP * W[I - 1]

        PQ[KPT - 1] = ZERO
        IP = trunc(PTSID[KPT - 1])
        IQ = trunc(NP * PTSID[KPT - 1] - (IP * NP))
        if (IP > 0):
            XP = PTSAUX0[IP - 1]
            XPT[KPT - 1][IP - 1] = XP

        if (IQ > 0):
            XQ = PTSAUX0[IQ - 1]
            if (IP == 0):
                XQ = PTSAUX1[IQ - 1]

            XPT[KPT - 1][IQ - 1] = XQ

        #
        #     Set VQUAD to the value of the current model at the new point.
        #
        VQUAD = FBASE
        if (IP > 0):
            IHP = trunc((IP + IP * IP) / 2)
            VQUAD = VQUAD + XP * (GOPT[IP - 1] + HALF * XP * HQ[IHP - 1])

        if (IQ > 0):
            IHQ = trunc((IQ + IQ * IQ) / 2)
            VQUAD = VQUAD + XQ * (GOPT[IQ - 1] + HALF * XQ * HQ[IHQ - 1])
            if (IP > 0):
                IW = trunc(max(IHP, IHQ)) - trunc(abs(IP - IQ))
                VQUAD = VQUAD + XP * XQ * HQ[IW - 1]

        for K in range(1, NPT+1):
            TEMP = ZERO
            if (IP > 0):
                TEMP = TEMP + XP * XPT[K - 1][IP - 1]

            if (IQ > 0):
                TEMP = TEMP + XQ * XPT[K - 1][IQ - 1]

            VQUAD = VQUAD + HALF * PQ[K - 1] * TEMP * TEMP

        #
        #     Calculate F at the new interpolation point, and set DIFF to the factor
        #     that is going to multiply the KPT-th Lagrange function when the model
        #     is updated to provide interpolation to the new function value.
        #
        for I in range(1, N+1):
            W[I - 1] = min(max(XL[I - 1], XBASE[I - 1] +
                           XPT[KPT - 1][I - 1]), XU[I - 1])
            if (XPT[KPT - 1][I - 1] == SL[I - 1]):
                W[I - 1] = XL[I - 1]

            if (XPT[KPT - 1][I - 1] == SU[I - 1]):
                W[I - 1] = XU[I - 1]

        NF = NF + 1
        F = calfun(N, W)
        if (IPRINT == 3):
            print(f'Function number: {NF}    F = {F}    ' +
                  f'The corresponding X is: {W[: N]}')

        FVAL[KPT - 1] = F
        if (F < FVAL[KOPT - 1]):
            KOPT = KPT

        DIFF = F - VQUAD
        #
        #     Update the quadratic model. The return from the subroutine occurs when
        #     all the new interpolation points are included in the model.
        #
        for I in range(1, N+1):
            GOPT[I - 1] = GOPT[I - 1] + DIFF * BMAT[KPT - 1][I - 1]

        for K in range(1, NPT+1):
            SUM = ZERO
            for J in range(1, NPTM+1):
                SUM = SUM + ZMAT[K - 1][J - 1] * ZMAT[KPT - 1][J - 1]

            TEMP = DIFF * SUM
            if (PTSID[K - 1] == ZERO):
                PQ[K - 1] = PQ[K - 1] + TEMP
            else:
                IP = trunc(PTSID[K - 1])
                IQ = trunc(NP * PTSID[K - 1] - (IP * NP))
                IHQ = trunc((IQ * IQ + IQ) / 2)
                if (IP == 0):
                    HQ[IHQ - 1] = HQ[IHQ - 1] + TEMP * \
                        PTSAUX1[IQ - 1] ** 2
                else:
                    IHP = trunc((IP * IP + IP) / 2)
                    HQ[IHP - 1] = HQ[IHP - 1] + TEMP * \
                        PTSAUX0[IP - 1] ** 2
                    if (IQ > 0):
                        HQ[IHQ - 1] = HQ[IHQ - 1] + TEMP * \
                            PTSAUX0[IQ - 1] ** 2
                        IW = trunc(max(IHP, IHQ)) - trunc(abs(IQ - IP))
                        HQ[IW - 1] = HQ[IW - 1] + TEMP * \
                            PTSAUX0[IP - 1] * PTSAUX0[IQ - 1]

        PTSID[KPT - 1] = ZERO

    return (XBASE, XPT, FVAL, XOPT, GOPT, HQ, PQ, BMAT, ZMAT, SL, SU, NF,
            KOPT, VLAG, PTSAUX0, PTSAUX1, PTSID, W)


def prelim(N, NPT, X, XL, XU, RHOBEG, IPRINT, MAXFUN, XBASE, XPT, FVAL,
           GOPT, HQ, PQ, BMAT, ZMAT, SL, SU):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*),
    #      1     HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
    #
    #     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
    #       same as the corresponding arguments in SUBROUTINE BOBYQA.
    #     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
    #       are the same as the corresponding arguments in BOBYQB, the elements
    #       of SL and SU being set in BOBYQA.
    #     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
    #       it is set by PRELIM to the gradient of the quadratic model at XBASE.
    #       If XOPT is nonzero, BOBYQB will change it to its usual value later.
    #     NF is maintaned as the number of calls of CALFUN so far.
    #     KOPT will be such that the least calculated value of F so far is at
    #       the point XPT(KOPT,.)+XBASE in the space of the variables.
    #
    #     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
    #     BMAT and ZMAT for the first iteration, and it maintains the values of
    #     NF and KOPT. The vector X is also changed by PRELIM.
    #
    #     Set some constants.
    #
    HALF = 0.5e0
    ONE = 1.0e0
    TWO = 2.0e0
    ZERO = 0.0e0
    RHOSQ = RHOBEG * RHOBEG
    RECIP = ONE / RHOSQ
    NP = N + 1
    #
    #     Set XBASE to the initial vector of variables, and set the initial
    #     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
    #     elements in matlab are all set to zeros() omit init
    #
    for J in range(1, N+1):
        XBASE[J - 1] = X[J - 1]

    #
    #     Begin the initialization procedure. NF becomes one more than the number
    #     of function values so far. The coordinates of the displacement of the
    #     next initial interpolation point from XBASE are set in XPT(NF+1,.).
    #
    NF = 0
    while (NF < NPT and NF < MAXFUN):
        NFM = NF
        NFX = NF - N
        NF = NF + 1
        if (NFM <= 2 * N):
            if (NFM >= 1 and NFM <= N):
                STEPA = RHOBEG
                if (SU[NFM - 1] == ZERO):
                    STEPA = -STEPA

                XPT[NF - 1][NFM - 1] = STEPA
            elif (NFM > N):
                STEPA = XPT[NF - N - 1][NFX - 1]
                STEPB = -RHOBEG
                if (SL[NFX - 1] == ZERO):
                    STEPB = min(TWO * RHOBEG, SU[NFX - 1])

                if (SU[NFX - 1] == ZERO):
                    STEPB = max(-TWO * RHOBEG, SL[NFX - 1])

                XPT[NF - 1][NFX - 1] = STEPB

        else:
            ITEMP = (NFM - NP) / N
            JPT = NFM - ITEMP * N - N
            IPT = JPT + ITEMP
            if (IPT > N):
                ITEMP = JPT
                JPT = IPT - N
                IPT = ITEMP

            XPT[NF - 1][IPT - 1] = XPT[IPT][IPT - 1]
            XPT[NF - 1][JPT - 1] = XPT[JPT][JPT - 1]

        #
        #     Calculate the next value of F. The least function value so far and
        #     its index are required.
        #
        for J in range(1, N+1):
            X[J - 1] = min(max(XL[J - 1], XBASE[J - 1] +
                           XPT[NF - 1][J - 1]), XU[J - 1])
            if (XPT[NF - 1][J - 1] == SL[J - 1]):
                X[J - 1] = XL[J - 1]

            if (XPT[NF - 1][J - 1] == SU[J - 1]):
                X[J - 1] = XU[J - 1]

        F = calfun(N, X)
        if (IPRINT == 3):
            print(f'Function number {NF}    F = {F}' +
                  f' The corresponding X is: {X[: N]}')

        FVAL[NF - 1] = F
        if (NF == 1):
            FBEG = F
            KOPT = 1
        elif (F < FVAL[KOPT - 1]):
            KOPT = NF

        #
        #     Set the nonzero initial elements of BMAT and the quadratic model in the
        #     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
        #     of the NF-th and (NF-N)-th interpolation points may be switched, in
        #     order that the function value at the first of them contributes to the
        #     off-diagonal second derivative terms of the initial quadratic model.
        #
        if (NF <= 2 * N + 1):
            if (NF >= 2 and NF <= N + 1):
                GOPT[NFM - 1] = (F - FBEG) / STEPA

                if (NPT < NF + N):
                    BMAT[0][NFM - 1] = -ONE / STEPA
                    BMAT[NF - 1][NFM - 1] = ONE / STEPA
                    BMAT[NPT + NFM - 1][NFM - 1] = -HALF * RHOSQ

            elif (NF >= N + 2):
                IH = trunc((NFX * (NFX + 1)) / 2)
                TEMP = (F - FBEG) / STEPB
                DIFF = STEPB - STEPA
                HQ[IH - 1] = TWO * (TEMP - GOPT[NFX - 1]) / DIFF
                GOPT[NFX - 1] = (GOPT[NFX - 1] * STEPB - TEMP * STEPA) / DIFF
                if (STEPA * STEPB < ZERO):
                    if (F < FVAL[NF - N - 1]):
                        FVAL[NF - 1] = FVAL[NF - N - 1]
                        FVAL[NF - N - 1] = F
                        if (KOPT == NF):
                            KOPT = NF - N

                        XPT[NF - N - 1][NFX - 1] = STEPB
                        XPT[NF - 1][NFX - 1] = STEPA

                BMAT[0][NFX - 1] = - (STEPA + STEPB) / (STEPA * STEPB)
                BMAT[NF - 1][NFX - 1] = -HALF / XPT[NF - N - 1][NFX - 1]
                BMAT[NF - N - 1][NFX - 1] = - \
                    BMAT[0][NFX - 1] - BMAT[NF - 1][NFX - 1]
                ZMAT[0][NFX - 1] = sqrt(TWO) / (STEPA * STEPB)
                ZMAT[NF - 1][NFX - 1] = sqrt(HALF) / RHOSQ
                ZMAT[NF - N - 1][NFX - 1] = - \
                    ZMAT[0][NFX - 1] - ZMAT[NF - 1][NFX - 1]

            #
            #     Set the off-diagonal second derivatives of the Lagrange functions and
            #     the initial quadratic model.
            #
        else:
            IH = (IPT * (IPT - 1)) / 2 + JPT
            ZMAT[0][NFX - 1] = RECIP
            ZMAT[NF - 1][NFX - 1] = RECIP
            ZMAT[IPT][NFX - 1] = -RECIP
            ZMAT[JPT][NFX - 1] = -RECIP
            TEMP = XPT[NF - 1][IPT - 1] * XPT[NF - 1][JPT - 1]
            HQ[IH - 1] = (FBEG - FVAL[IPT] - FVAL[JPT] + F) / TEMP

    return (X, XBASE, XPT, FVAL, GOPT, HQ, PQ, BMAT, ZMAT, NF, KOPT)


def calfun(N, X):
    F = 0.0e0
    for I in range(4, N+1, 2):
        for J in range(2, (I - 1), 2):
            TEMP = (X[I - 2] - X[J - 2]) ** 2 + (X[I-1] - X[J-1]) ** 2
            TEMP = max(TEMP, 1.0e-6)
            F = F + 1 / sqrt(TEMP)
    return F


def altmov(N, NPT, XPT, XOPT, BMAT, ZMAT, SL, SU, KOPT, KNEW, ADELT, XNEW, XALT, CAUCHY, GLAG, HCOL, W):
    #          IMPLICIT REAL*8 (A-H,O-Z)
    #          DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),
    #      1     SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
    #
    #     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
    #       the same meanings as the corresponding arguments of BOBYQB.
    #     KOPT is the index of the optimal interpolation point.
    #     KNEW is the index of the interpolation point that is going to be moved.
    #     ADELT is the current trust region bound.
    #     XNEW will be set to a suitable new position for the interpolation point
    #       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
    #       bounds and it should provide a large denominator in the next call of
    #       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
    #       straight lines through XOPT and another interpolation point.
    #     XALT also provides a large value of the modulus of the KNEW-th Lagrange
    #       function subject to the constraints that have been mentioned, its main
    #       difference from XNEW being that XALT-XOPT is a constrained version of
    #       the Cauchy step within the trust region. An exception is that XALT is
    #       not calculated if all components of GLAG (see below) are zero.
    #     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
    #     CAUCHY will be set to the square of the KNEW-th Lagrange function at
    #       the step XALT-XOPT from XOPT for the vector XALT that is returned,
    #       except that CAUCHY is set to zero if XALT is not calculated.
    #     GLAG is a working space vector of length N for the gradient of the
    #       KNEW-th Lagrange function at XOPT.
    #     HCOL is a working space vector of length NPT for the second derivative
    #       coefficients of the KNEW-th Lagrange function.
    #     W is a working space vector of length 2N that is going to hold the
    #       constrained Cauchy step from XOPT of the Lagrange function, followed
    #       by the downhill version of XALT when the uphill step is calculated.
    #
    #     Set the first NPT components of W to the leading elements of the
    #     KNEW-th column of the H matrix.
    #
    HALF = 0.5e0
    ONE = 1.0e0
    ZERO = 0.0e0
    CONST = ONE + sqrt(2.0e0)
    for K in range(1, NPT+1):
        HCOL[K - 1] = ZERO

    for J in range(1, NPT - N):
        TEMP = ZMAT[KNEW - 1][J - 1]
        for K in range(1, NPT+1):
            HCOL[K - 1] = HCOL[K - 1] + TEMP * ZMAT[K - 1][J - 1]

    ALPHA = HCOL[KNEW - 1]
    HA = HALF * ALPHA
    #
    #     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
    #
    for I in range(1, N+1):
        GLAG[I - 1] = BMAT[KNEW - 1][I - 1]

    for K in range(1, NPT+1):
        TEMP = ZERO
        for J in range(1, N+1):
            TEMP = TEMP + XPT[K - 1][J - 1] * XOPT[J - 1]

        TEMP = HCOL[K - 1] * TEMP
        for I in range(1, N+1):
            GLAG[I - 1] = GLAG[I - 1] + TEMP * XPT[K - 1][I - 1]

    #
    #     Search for a large denominator along the straight lines through XOPT
    #     and another interpolation point. SLBD and SUBD will be lower and upper
    #     bounds on the step along each of these lines in turn. PREDSQ will be
    #     set to the square of the predicted denominator for each line. PRESAV
    #     will be set to the largest admissible value of PREDSQ that occurs.
    #
    PRESAV = ZERO
    for K in range(1, NPT+1):
        if (K == KOPT):
            continue

        DDERIV = ZERO
        DISTSQ = ZERO
        for I in range(1, N+1):
            TEMP = XPT[K - 1][I - 1] - XOPT[I - 1]
            DDERIV = DDERIV + GLAG[I - 1] * TEMP
            DISTSQ = DISTSQ + TEMP * TEMP

        SUBD = ADELT / sqrt(DISTSQ)
        SLBD = -SUBD
        ILBD = 0
        IUBD = 0
        SUMIN = min(ONE, SUBD)
        #
        #     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
        #
        for I in range(1, N+1):
            TEMP = XPT[K - 1][I - 1] - XOPT[I - 1]
            if (TEMP > ZERO):
                if (SLBD * TEMP < SL[I - 1] - XOPT[I - 1]):
                    SLBD = (SL[I - 1] - XOPT[I - 1]) / TEMP
                    ILBD = -I

                if (SUBD * TEMP > SU[I - 1] - XOPT[I - 1]):
                    SUBD = max(SUMIN, (SU[I - 1] - XOPT[I - 1]) / TEMP)
                    IUBD = I

            elif (TEMP < ZERO):
                if (SLBD * TEMP > SU[I - 1] - XOPT[I - 1]):
                    SLBD = (SU[I - 1] - XOPT[I - 1]) / TEMP
                    ILBD = I

                if (SUBD * TEMP < SL[I - 1] - XOPT[I - 1]):
                    SUBD = max(SUMIN, (SL[I - 1] - XOPT[I - 1]) / TEMP)
                    IUBD = -I

        #
        #     Seek a large modulus of the KNEW-th Lagrange function when the index
        #     of the other interpolation point on the line through XOPT is KNEW.
        #
        if (K == KNEW):
            DIFF = DDERIV - ONE
            STEP = SLBD
            VLAG = SLBD * (DDERIV - SLBD * DIFF)
            ISBD = ILBD
            TEMP = SUBD * (DDERIV - SUBD * DIFF)
            if (abs(TEMP) > abs(VLAG)):
                STEP = SUBD
                VLAG = TEMP
                ISBD = IUBD

            TEMPD = HALF * DDERIV
            TEMPA = TEMPD - DIFF * SLBD
            TEMPB = TEMPD - DIFF * SUBD
            if (TEMPA * TEMPB < ZERO):
                TEMP = TEMPD * TEMPD / DIFF
                if (abs(TEMP) > abs(VLAG)):
                    STEP = TEMPD / DIFF
                    VLAG = TEMP
                    ISBD = 0

            #
            #     Search along each of the other lines through XOPT and another point.
            #
        else:
            STEP = SLBD
            VLAG = SLBD * (ONE - SLBD)
            ISBD = ILBD
            TEMP = SUBD * (ONE - SUBD)
            if (abs(TEMP) > abs(VLAG)):
                STEP = SUBD
                VLAG = TEMP
                ISBD = IUBD

            if (SUBD > HALF):
                if (abs(VLAG) < 0.25e0):
                    STEP = HALF
                    VLAG = 0.25e0
                    ISBD = 0

            VLAG = VLAG * DDERIV

        #
        #     Calculate PREDSQ for the current line search and maintain PRESAV.
        #
        TEMP = STEP * (ONE - STEP) * DISTSQ
        PREDSQ = VLAG * VLAG * (VLAG * VLAG + HA * TEMP * TEMP)
        if (PREDSQ > PRESAV):
            PRESAV = PREDSQ
            KSAV = K
            STPSAV = STEP
            IBDSAV = ISBD

    #
    #     Construct XNEW in a way that satisfies the bound constraints exactly.
    #
    for I in range(1, N+1):
        TEMP = XOPT[I - 1] + STPSAV * (XPT[KSAV - 1][I - 1] - XOPT[I - 1])
        XNEW[I - 1] = max(SL[I - 1], min(SU[I - 1], TEMP))

    if (IBDSAV < 0):
        XNEW[-IBDSAV - 1] = SL[-IBDSAV - 1]

    if (IBDSAV > 0):
        XNEW[IBDSAV - 1] = SU[IBDSAV - 1]

    #
    #     Prepare for the iterative method that assembles the constrained Cauchy
    #     step in W. The sum of squares of the fixed components of W is formed in
    #     WFIXSQ, and the free components of W are set to BIGSTP.
    #
    BIGSTP = ADELT + ADELT
    IFLAG = 0
    while (True):
        WFIXSQ = ZERO
        GGFREE = ZERO
        for I in range(1, N+1):
            W[I - 1] = ZERO
            TEMPA = min(XOPT[I - 1] - SL[I - 1], GLAG[I - 1])
            TEMPB = max(XOPT[I - 1] - SU[I - 1], GLAG[I - 1])
            if (TEMPA > ZERO or TEMPB < ZERO):
                W[I - 1] = BIGSTP
                GGFREE = GGFREE + GLAG[I - 1] ** 2

        if (GGFREE == ZERO):
            CAUCHY = ZERO
            return (XNEW, XALT, ALPHA, CAUCHY, GLAG, HCOL, W)

        #
        #     Investigate whether more components of W can be fixed.
        #
        while (True):
            TEMP = ADELT * ADELT - WFIXSQ
            if (TEMP > ZERO):
                WSQSAV = WFIXSQ
                STEP = sqrt(TEMP / GGFREE)
                GGFREE = ZERO
                for I in range(1, N+1):
                    if (W[I - 1] == BIGSTP):
                        TEMP = XOPT[I - 1] - STEP * GLAG[I - 1]
                        if (TEMP <= SL[I - 1]):
                            W[I - 1] = SL[I - 1] - XOPT[I - 1]
                            WFIXSQ = WFIXSQ + W[I - 1] ** 2
                        elif (TEMP >= SU[I - 1]):
                            W[I - 1] = SU[I - 1] - XOPT[I - 1]
                            WFIXSQ = WFIXSQ + W[I - 1] ** 2
                        else:
                            GGFREE = GGFREE + GLAG[I - 1] ** 2

                if (WFIXSQ <= WSQSAV or GGFREE <= ZERO):
                    break

            else:
                break

        #
        #     Set the remaining free components of W and all components of XALT,
        #     except that W may be scaled later.
        #
        GW = ZERO
        for I in range(1, N+1):
            if (W[I - 1] == BIGSTP):
                W[I - 1] = -STEP * GLAG[I - 1]
                XALT[I - 1] = max(SL[I - 1], min(SU[I - 1],
                                  XOPT[I - 1] + W[I - 1]))
            elif (W[I - 1] == ZERO):
                XALT[I - 1] = XOPT[I - 1]
            elif (GLAG[I - 1] > ZERO):
                XALT[I - 1] = SL[I - 1]
            else:
                XALT[I - 1] = SU[I - 1]

            GW = GW + GLAG[I - 1] * W[I - 1]

        #
        #     Set CURV to the curvature of the KNEW-th Lagrange function along W.
        #     Scale W by a factor less than one if that can reduce the modulus of
        #     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
        #     the square of this function.
        #
        CURV = ZERO
        for K in range(1, NPT+1):
            TEMP = ZERO
            for J in range(1, N+1):
                TEMP = TEMP + XPT[K - 1][J - 1] * W[J - 1]

            CURV = CURV + HCOL[K - 1] * TEMP * TEMP

        if (IFLAG == 1):
            CURV = -CURV

        if (CURV > -GW and CURV < -CONST * GW):
            SCALE = -GW / CURV
            for I in range(1, N+1):
                TEMP = XOPT[I - 1] + SCALE * W[I - 1]
                XALT[I - 1] = max(SL[I - 1], min(SU[I - 1], TEMP))

            CAUCHY = (HALF * GW * SCALE) ** 2
        else:
            CAUCHY = (GW + HALF * CURV) ** 2

        #
        #     If IFLAG is zero, then XALT is calculated as before after reversing
        #     the sign of GLAG. Thus two XALT vectors become available. The one that
        #     is chosen is the one that gives the larger value of CAUCHY.
        #
        if (IFLAG != 0):
            break

        for I in range(1, N+1):
            GLAG[I - 1] = -GLAG[I - 1]
            W[N + I - 1] = XALT[I - 1]

        CSAVE = CAUCHY
        IFLAG = 1

    if (CSAVE > CAUCHY):
        for I in range(1, N+1):
            XALT[I - 1] = W[N + I - 1]

        CAUCHY = CSAVE

    return (XNEW, XALT, ALPHA, CAUCHY, GLAG, HCOL, W)


if __name__ == '__main__':
    #
    #     Test problem for BOBYQA, the objective function being the sum of
    #     the reciprocals of all pairwise distances between the points P_I,
    #     I=1,2,...,M in two dimensions, where M=N/2 and where the components
    #     of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
    #     defines the M points P_I. The initial X gives equally spaced points
    #     on a circle. Four different choices of the pairs (N,NPT) are tried,
    #     namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
    #     minimum that is not global occurs in both the N=10 cases. The details
    #     of the results are highly sensitive to computer rounding errors. The
    #     choice IPRINT=2 provides the current X and optimal F so far whenever
    #     RHO is reduced. The bound constraints of the problem require every
    #     component of X to be in the interval [-1,1].
    #
    # IMPLICIT REAL * 8 (A - H, O - Z)
    # DIMENSION X(100), XL(100), XU(100), W(500000)
    #
    X = [0 for _ in range(100)]
    XL = [0 for _ in range(100)]
    XU = [0 for _ in range(100)]
    TWOPI = 8.0e0 * atan(1.0e0)
    BDL = -1.0e0
    BDU = 1.0e0
    IPRINT = 2
    MAXFUN = 500000
    RHOBEG = 1.0e-1
    RHOEND = 1.0e-6
    M = 10
    while (M <= 10):
        N = 2 * M
        for I in range(1, N+1):
            XL[I - 1] = BDL
            XU[I - 1] = BDU

        for JCASE in range(1, 2):
            NPT = N + 6
            if (JCASE == 2):
                NPT = 2 * N + 1

            print(f'\n[new]:2D output with M = {M}' +
                  f',  N = {N}' + f'  and  NPT = {NPT}.')
            for J in range(1, M+1):
                TEMP = J * TWOPI / M
                X[2 * J - 2] = cos(TEMP)
                X[2 * J - 1] = sin(TEMP)

            bobyqa(N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT, MAXFUN)

        M = M + M
