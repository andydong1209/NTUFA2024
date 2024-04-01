Global Const Pi = 3.14159265358979

Option Explicit     'Requirs that all variables to be declared explicitly.
Option Compare Text 'Uppercase letters to be equivalent to lowercase letters.

Option Base 1       'The "Option Base" statment alowws to specify 0 or 1 as the
                    'default first index of arrays.
                    


Public Function Max(X, y)
            Max = Application.Max(X, y)
End Function

Public Function Min(X, y)
            Min = Application.Min(X, y)
End Function
                                                                                                         
'// The normal distribution function
Public Function ND(X As Double) As Double
    ND = 1 / Sqr(2 * Pi) * Exp(-X ^ 2 / 2)
End Function


'// The cumulative normal distribution function
Public Function CND(X As Double) As Double
    
    Dim L As Double, K As Double
    Const a1 = 0.31938153:  Const a2 = -0.356563782: Const a3 = 1.781477937:
    Const a4 = -1.821255978:  Const a5 = 1.330274429
    
    L = Abs(X)
    K = 1 / (1 + 0.2316419 * L)
    CND = 1 - 1 / Sqr(2 * Pi) * Exp(-L ^ 2 / 2) * (a1 * K + a2 * K ^ 2 + a3 * K ^ 3 + a4 * K ^ 4 + a5 * K ^ 5)
    
    If X < 0 Then
        CND = 1 - CND
    End If
End Function


'// The cumulative bivariate normal distribution function
Public Function CBND(a As Double, b As Double, rho As Double) As Double

    Dim X As Variant, y As Variant
    Dim rho1 As Double, rho2 As Double, delta As Double
    Dim a1 As Double, b1 As Double, Sum As Double
    Dim I As Integer, j As Integer
    
    X = Array(0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334)
    y = Array(0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604)
    a1 = a / Sqr(2 * (1 - rho ^ 2))
    b1 = b / Sqr(2 * (1 - rho ^ 2))
    
    If a <= 0 And b <= 0 And rho <= 0 Then
        Sum = 0
        For I = 1 To 5
            For j = 1 To 5
                Sum = Sum + X(I) * X(j) * Exp(a1 * (2 * y(I) - a1) _
                + b1 * (2 * y(j) - b1) + 2 * rho * (y(I) - a1) * (y(j) - b1))
            Next
        Next
        CBND = Sqr(1 - rho ^ 2) / Pi * Sum
    ElseIf a <= 0 And b >= 0 And rho >= 0 Then
        CBND = CND(a) - CBND(a, -b, -rho)
    ElseIf a >= 0 And b <= 0 And rho >= 0 Then
        CBND = CND(b) - CBND(-a, b, -rho)
    ElseIf a >= 0 And b >= 0 And rho <= 0 Then
        CBND = CND(a) + CND(b) - 1 + CBND(-a, -b, rho)
    ElseIf a * b * rho > 0 Then
        rho1 = (rho * a - b) * Sgn(a) / Sqr(a ^ 2 - 2 * rho * a * b + b ^ 2)
        rho2 = (rho * b - a) * Sgn(b) / Sqr(a ^ 2 - 2 * rho * a * b + b ^ 2)
        delta = (1 - Sgn(a) * Sgn(b)) / 4
        CBND = CBND(a, 0, rho1) + CBND(b, 0, rho2) - delta
    End If
End Function


'// Black and Scholes (1973) Stock options
Public Function BlackScholes(CallPutFlag As String, S As Double, X _
                As Double, T As Double, r As Double, v As Double) As Double
                
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (r + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    If CallPutFlag = "c" Then
        BlackScholes = S * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        BlackScholes = X * Exp(-r * T) * CND(-d2) - S * CND(-d1)
    End If
End Function


'// Merton (1973) Options on stock indices
Public Function Merton73(CallPutFlag As String, S As Double, X _
                As Double, T As Double, r As Double, q As Double, v As Double) As Double
                
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (r - q + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    If CallPutFlag = "c" Then
        Merton73 = S * Exp(-q * T) * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        Merton73 = X * Exp(-r * T) * CND(-d2) - S * Exp(-q * T) * CND(-d1)
    End If
End Function


'// Black (1977) Options on futures/forwards
Public Function Black76(CallPutFlag As String, F As Double, X _
                As Double, T As Double, r As Double, v As Double) As Double
                
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(F / X) + (v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    If CallPutFlag = "c" Then
        Black76 = Exp(-r * T) * (F * CND(d1) - X * CND(d2))
    ElseIf CallPutFlag = "p" Then
        Black76 = Exp(-r * T) * (X * CND(-d2) - F * CND(-d1))
    End If
End Function


'// Garman and Kohlhagen (1983) Currency options
Public Function GarmanKolhagen(CallPutFlag As String, S As Double, X _
                As Double, T As Double, r As Double, rf As Double, v As Double) As Double
                
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (r - rf + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    If CallPutFlag = "c" Then
        GarmanKolhagen = S * Exp(-rf * T) * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        GarmanKolhagen = X * Exp(-r * T) * CND(-d2) - S * Exp(-rf * T) * CND(-d1)
    End If
End Function


'// The generalized Black and Scholes formula
Public Function GBlackScholes(CallPutFlag As String, S As Double, X _
                As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)

    If CallPutFlag = "c" Then
        GBlackScholes = S * Exp((b - r) * T) * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        GBlackScholes = X * Exp(-r * T) * CND(-d2) - S * Exp((b - r) * T) * CND(-d1)
    End If
End Function


'// Delta for the generalized Black and Scholes formula
Public Function GDelta(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, _
                b As Double, v As Double) As Double
                
    Dim d1 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    
    If CallPutFlag = "c" Then
        GDelta = Exp((b - r) * T) * CND(d1)
    ElseIf CallPutFlag = "p" Then
        GDelta = Exp((b - r) * T) * (CND(d1) - 1)
    End If
End Function


'// Gamma for the generalized Black and Scholes formula
Public Function GGamma(S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    GGamma = Exp((b - r) * T) * ND(d1) / (S * v * Sqr(T))
End Function


'// Theta for the generalized Black and Scholes formula
Public Function GTheta(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)

    If CallPutFlag = "c" Then
        GTheta = -S * Exp((b - r) * T) * ND(d1) * v / (2 * Sqr(T)) - (b - r) * S * Exp((b - r) * T) * CND(d1) - r * X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        GTheta = -S * Exp((b - r) * T) * ND(d1) * v / (2 * Sqr(T)) + (b - r) * S * Exp((b - r) * T) * CND(-d1) + r * X * Exp(-r * T) * CND(-d2)
    End If
End Function


'// Vega for the generalized Black and Scholes formula
Public Function GVega(S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    GVega = S * Exp((b - r) * T) * ND(d1) * Sqr(T)
End Function


'// Rho for the generalized Black and Scholes formula
Public Function GRho(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
   Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    If CallPutFlag = "c" Then
        If b <> 0 Then
            GRho = T * X * Exp(-r * T) * CND(d2)
        Else
            GRho = -T * GBlackScholes(CallPutFlag, S, X, T, r, b, v)
        End If
    ElseIf CallPutFlag = "p" Then
        If b <> 0 Then
            GRho = -T * X * Exp(-r * T) * CND(-d2)
        Else
            GRho = -T * GBlackScholes(CallPutFlag, S, X, T, r, b, v)
        End If
    End If
End Function

'// Carry for the generalized Black and Scholes formula
Public Function GCarry(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double

    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    If CallPutFlag = "c" Then
        GCarry = T * S * Exp((b - r) * T) * CND(d1)
    ElseIf CallPutFlag = "p" Then
        GCarry = -T * S * Exp((b - r) * T) * CND(-d1)
    End If
End Function


'// French (1984) adjusted Black and Scholes model for trading day volatility
Public Function French(CallPutFlag As String, S As Double, X As Double, T As Double, t1 As Double, _
                r As Double, b As Double, v As Double) As Double

    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / X) + b * T + v ^ 2 / 2 * t1) / (v * Sqr(t1))
    d2 = d1 - v * Sqr(t1)
  
    If CallPutFlag = "c" Then
        French = S * Exp((b - r) * T) * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        French = X * Exp(-r * T) * CND(-d2) - S * Exp((b - r) * T) * CND(-d1)
    End If
End Function


'// Merton (1976) jump diffusion model
Public Function JumpDiffusion(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, v As Double, _
                lambda As Double, gamma As Double) As Double

    Dim delta As Double, Sum As Double
    Dim Z As Double, vi As Double
    Dim I As Integer

    delta = Sqr(gamma * v ^ 2 / lambda)
    Z = Sqr(v ^ 2 - lambda * delta ^ 2)
    Sum = 0
    For I = 0 To 10
        vi = Sqr(Z ^ 2 + delta ^ 2 * (I / T))
        Sum = Sum + Exp(-lambda * T) * (lambda * T) ^ I / Application.Fact(I) * _
        GBlackScholes(CallPutFlag, S, X, T, r, r, vi)
    Next
        JumpDiffusion = Sum
End Function


'// Miltersen Schwartz (1997) commodity option model
Public Function MiltersenSwartz(CallPutFlag As String, Pt As Double, FT As Double, X As Double, t1 As Double, _
                T2 As Double, vS As Double, vE As Double, vf As Double, rhoSe As Double, _
                rhoSf As Double, rhoef As Double, Kappae As Double, Kappaf As Double) As Double
    
                Dim vz As Double, vxz As Double
                Dim d1 As Double, d2 As Double
                
                vz = vS ^ 2 * t1 + 2 * vS * (vf * rhoSf * 1 / Kappaf * (t1 - 1 / Kappaf * Exp(-Kappaf * T2) * (Exp(Kappaf * t1) - 1)) _
                    - vE * rhoSe * 1 / Kappae * (t1 - 1 / Kappae * Exp(-Kappae * T2) * (Exp(Kappae * t1) - 1))) _
                    + vE ^ 2 * 1 / Kappae ^ 2 * (t1 + 1 / (2 * Kappae) * Exp(-2 * Kappae * T2) * (Exp(2 * Kappae * t1) - 1) - 2 * 1 / Kappae * Exp(-Kappae * T2) * (Exp(Kappae * t1) - 1)) _
                    + vf ^ 2 * 1 / Kappaf ^ 2 * (t1 + 1 / (2 * Kappaf) * Exp(-2 * Kappaf * T2) * (Exp(2 * Kappaf * t1) - 1) - 2 * 1 / Kappaf * Exp(-Kappaf * T2) * (Exp(Kappaf * t1) - 1)) _
                    - 2 * vE * vf * rhoef * 1 / Kappae * 1 / Kappaf * (t1 - 1 / Kappae * Exp(-Kappae * T2) * (Exp(Kappae * t1) - 1) - 1 / Kappaf * Exp(-Kappaf * T2) * (Exp(Kappaf * t1) - 1) _
                    + 1 / (Kappae + Kappaf) * Exp(-(Kappae + Kappaf) * T2) * (Exp((Kappae + Kappaf) * t1) - 1))
                
                vxz = vf * 1 / Kappaf * (vS * rhoSf * (t1 - 1 / Kappaf * (1 - Exp(-Kappaf * t1))) _
                    + vf * 1 / Kappaf * (t1 - 1 / Kappaf * Exp(-Kappaf * T2) * (Exp(Kappaf * t1) - 1) - 1 / Kappaf * (1 - Exp(-Kappaf * t1)) _
                    + 1 / (2 * Kappaf) * Exp(-Kappaf * T2) * (Exp(Kappaf * t1) - Exp(-Kappaf * t1))) _
                    - vE * rhoef * 1 / Kappae * (t1 - 1 / Kappae * Exp(-Kappae * T2) * (Exp(Kappae * t1) - 1) - 1 / Kappaf * (1 - Exp(-Kappaf * t1)) _
                    + 1 / (Kappae + Kappaf) * Exp(-Kappae * T2) * (Exp(Kappae * t1) - Exp(-Kappaf * t1))))
                
                vz = Sqr(vz)
                
                d1 = (Log(FT / X) - vxz + vz ^ 2 / 2) / vz
                d2 = (Log(FT / X) - vxz - vz ^ 2 / 2) / vz
                
                If CallPutFlag = "c" Then
                    MiltersenSwartz = Pt * (FT * Exp(-vxz) * CND(d1) - X * CND(d2))
                ElseIf CallPutFlag = "p" Then
                    MiltersenSwartz = Pt * (X * CND(-d2) - FT * Exp(-vxz) * CND(-d1))
                End If
                
End Function


'// * American options *


'// American Calls on stocks with known dividends, Roll-Geske-Whaley
Public Function RollGeskeWhaley(S As Double, X As Double, t1 As Double, T2 As Double, r As Double, d As Double, v As Double) As Double
    't1 time to dividend payout
    'T2 time to option expiration
    
    Dim Sx As Double, I As Double
    Dim a1 As Double, a2 As Double, b1 As Double, b2 As Double
    Dim HighS As Double, LowS As Double, epsilon As Double
    Dim ci As Double, infinity As Double
    
    infinity = 100000000
    epsilon = 0.00001
    Sx = S - d * Exp(-r * t1)
    If d <= X * (1 - Exp(-r * (T2 - t1))) Then '// Not optimal to exercise
        RollGeskeWhaley = BlackScholes("c", Sx, X, T2, r, v)
        Exit Function
    End If
    ci = BlackScholes("c", S, X, T2 - t1, r, v)
    HighS = S
    While (ci - HighS - d + X) > 0 And HighS < infinity
        HighS = HighS * 2
        ci = BlackScholes("c", HighS, X, T2 - t1, r, v)
    Wend
    If HighS > infinity Then
        RollGeskeWhaley = BlackScholes("c", Sx, X, T2, r, v)
        Exit Function
    End If
    
    LowS = 0
    I = HighS * 0.5
    ci = BlackScholes("c", I, X, T2 - t1, r, v)
    
    '// Search algorithm to find the critical stock price I
    While Abs(ci - I - d + X) > epsilon And HighS - LowS > epsilon
        If (ci - I - d + X) < 0 Then
            HighS = I
        Else
            LowS = I
        End If
        I = (HighS + LowS) / 2
        ci = BlackScholes("c", I, X, T2 - t1, r, v)
    Wend
    a1 = (Log(Sx / X) + (r + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    a2 = a1 - v * Sqr(T2)
    b1 = (Log(Sx / I) + (r + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    b2 = b1 - v * Sqr(t1)
   
    RollGeskeWhaley = Sx * CND(b1) + Sx * CBND(a1, -b1, -Sqr(t1 / T2)) - X * Exp(-r * T2) * CBND(a2, -b2, -Sqr(t1 / T2)) - (X - d) * Exp(-r * t1) * CND(b2)
End Function


'// The Barone-Adesi and Whaley (1987) American approximation
Public Function BAWAmericanApprox(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    If CallPutFlag = "c" Then
        BAWAmericanApprox = BAWAmericanCallApprox(S, X, T, r, b, v)
    ElseIf CallPutFlag = "p" Then
        BAWAmericanApprox = BAWAmericanPutApprox(S, X, T, r, b, v)
    End If
End Function

'// American call
Private Function BAWAmericanCallApprox(S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim Sk As Double, n As Double, K As Double
    Dim d1 As Double, Q2 As Double, a2 As Double

    If b >= r Then
        BAWAmericanCallApprox = GBlackScholes("c", S, X, T, r, b, v)
    Else
        Sk = Kc(X, T, r, b, v)
        n = 2 * b / v ^ 2                                           '
        K = 2 * r / (v ^ 2 * (1 - Exp(-r * T)))
        d1 = (Log(Sk / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
        Q2 = (-(n - 1) + Sqr((n - 1) ^ 2 + 4 * K)) / 2
        a2 = (Sk / Q2) * (1 - Exp((b - r) * T) * CND(d1))
        If S < Sk Then
            BAWAmericanCallApprox = GBlackScholes("c", S, X, T, r, b, v) + a2 * (S / Sk) ^ Q2
        Else
            BAWAmericanCallApprox = S - X
        End If
    End If
End Function

'// Newton Raphson algorithm to solve for the critical commodity price for a Call
Private Function Kc(X As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim n As Double, m As Double
    Dim Su As Double, Si As Double
    Dim h2 As Double, K As Double
    Dim d1 As Double, Q2 As Double, q2u As Double
    Dim LHS As Double, RHS As Double
    Dim bi As Double, E As Double
    
    '// Calculation of seed value, Si
    n = 2 * b / v ^ 2
    m = 2 * r / v ^ 2
    q2u = (-(n - 1) + Sqr((n - 1) ^ 2 + 4 * m)) / 2
    Su = X / (1 - 1 / q2u)
    h2 = -(b * T + 2 * v * Sqr(T)) * X / (Su - X)
    Si = X + (Su - X) * (1 - Exp(h2))

    K = 2 * r / (v ^ 2 * (1 - Exp(-r * T)))
    d1 = (Log(Si / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    Q2 = (-(n - 1) + Sqr((n - 1) ^ 2 + 4 * K)) / 2
    LHS = Si - X
    RHS = GBlackScholes("c", Si, X, T, r, b, v) + (1 - Exp((b - r) * T) * CND(d1)) * Si / Q2
    bi = Exp((b - r) * T) * CND(d1) * (1 - 1 / Q2) + (1 - Exp((b - r) * T) * CND(d1) / (v * Sqr(T))) / Q2
    E = 0.000001
    '// Newton Raphson algorithm for finding critical price Si
    While Abs(LHS - RHS) / X > E
        Si = (X + RHS - bi * Si) / (1 - bi)
        d1 = (Log(Si / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
        LHS = Si - X
        RHS = GBlackScholes("c", Si, X, T, r, b, v) + (1 - Exp((b - r) * T) * CND(d1)) * Si / Q2
        bi = Exp((b - r) * T) * CND(d1) * (1 - 1 / Q2) + (1 - Exp((b - r) * T) * ND(d1) / (v * Sqr(T))) / Q2
    Wend
        Kc = Si
End Function
'// American put
Private Function BAWAmericanPutApprox(S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim Sk As Double, n As Double, K As Double
    Dim d1 As Double, Q1 As Double, a1 As Double

    Sk = Kp(X, T, r, b, v)
    n = 2 * b / v ^ 2
    K = 2 * r / (v ^ 2 * (1 - Exp(-r * T)))
    d1 = (Log(Sk / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    Q1 = (-(n - 1) - Sqr((n - 1) ^ 2 + 4 * K)) / 2
    a1 = -(Sk / Q1) * (1 - Exp((b - r) * T) * CND(-d1))

    If S > Sk Then
        BAWAmericanPutApprox = GBlackScholes("p", S, X, T, r, b, v) + a1 * (S / Sk) ^ Q1
    Else
        BAWAmericanPutApprox = X - S
    End If
End Function

'// Newton Raphson algorithm to solve for the critical commodity price for a Put
Private Function Kp(X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
   
    Dim n As Double, m As Double
    Dim Su As Double, Si As Double
    Dim h1 As Double, K As Double
    Dim d1 As Double, q1u As Double, Q1 As Double
    Dim LHS As Double, RHS As Double
    Dim bi As Double, E As Double
    
    '// Calculation of seed value, Si
    n = 2 * b / v ^ 2
    m = 2 * r / v ^ 2
    q1u = (-(n - 1) - Sqr((n - 1) ^ 2 + 4 * m)) / 2
    Su = X / (1 - 1 / q1u)
    h1 = (b * T - 2 * v * Sqr(T)) * X / (X - Su)
    Si = Su + (X - Su) * Exp(h1)

    
    K = 2 * r / (v ^ 2 * (1 - Exp(-r * T)))
    d1 = (Log(Si / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    Q1 = (-(n - 1) - Sqr((n - 1) ^ 2 + 4 * K)) / 2
    LHS = X - Si
    RHS = GBlackScholes("p", Si, X, T, r, b, v) - (1 - Exp((b - r) * T) * CND(-d1)) * Si / Q1
    bi = -Exp((b - r) * T) * CND(-d1) * (1 - 1 / Q1) - (1 + Exp((b - r) * T) * ND(-d1) / (v * Sqr(T))) / Q1
    E = 0.000001
    '// Newton Raphson algorithm for finding critical price Si
    While Abs(LHS - RHS) / X > E
        Si = (X - RHS + bi * Si) / (1 + bi)
        d1 = (Log(Si / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
        LHS = X - Si
        RHS = GBlackScholes("p", Si, X, T, r, b, v) - (1 - Exp((b - r) * T) * CND(-d1)) * Si / Q1
        bi = -Exp((b - r) * T) * CND(-d1) * (1 - 1 / Q1) - (1 + Exp((b - r) * T) * CND(-d1) / (v * Sqr(T))) / Q1
    Wend
    Kp = Si
End Function


'// The Bjerksund and Stensland (1993) American approximation
Public Function BSAmericanApprox(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    If CallPutFlag = "c" Then
        BSAmericanApprox = BSAmericanCallApprox(S, X, T, r, b, v)
    ElseIf CallPutFlag = "p" Then  '// Use the Bjerksund and Stensland put-call transformation
        BSAmericanApprox = BSAmericanCallApprox(X, S, T, r - b, -b, v)
    End If
End Function

Public Function BSAmericanCallApprox(S As Double, X As Double, T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim BInfinity As Double, B0 As Double
    Dim ht As Double, I As Double
    Dim alpha As Double, Beta As Double
    
    If b >= r Then '// Never optimal to exersice before maturity
            BSAmericanCallApprox = GBlackScholes("c", S, X, T, r, b, v)
    Else
        Beta = (1 / 2 - b / v ^ 2) + Sqr((b / v ^ 2 - 1 / 2) ^ 2 + 2 * r / v ^ 2)
        BInfinity = Beta / (Beta - 1) * X
        B0 = Max(X, r / (r - b) * X)
        ht = -(b * T + 2 * v * Sqr(T)) * B0 / (BInfinity - B0)
        I = B0 + (BInfinity - B0) * (1 - Exp(ht))
        alpha = (I - X) * I ^ (-Beta)
        If S >= I Then
            BSAmericanCallApprox = S - X
        Else
            BSAmericanCallApprox = alpha * S ^ Beta - alpha * phi(S, T, Beta, I, I, r, b, v) + phi(S, T, 1, I, I, r, b, v) - phi(S, T, 1, X, I, r, b, v) - X * phi(S, T, 0, I, I, r, b, v) + X * phi(S, T, 0, X, I, r, b, v)
        End If
    End If
End Function

Private Function phi(S As Double, T As Double, gamma As Double, H As Double, I As Double, _
        r As Double, b As Double, v As Double) As Double

    Dim lambda As Double, kappa As Double
    Dim d As Double
    
    lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1) * v ^ 2) * T
    d = -(Log(S / H) + (b + (gamma - 0.5) * v ^ 2) * T) / (v * Sqr(T))
    kappa = 2 * b / (v ^ 2) + (2 * gamma - 1)
    phi = Exp(lambda) * S ^ gamma * (CND(d) - (I / S) ^ kappa * CND(d - 2 * Log(I / S) / (v * Sqr(T))))
End Function


'// Executive stock options
Public Function Executive(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, b As Double, v As Double, lambda As Double) As Double
    
    Executive = Exp(-lambda * T) * GBlackScholes(CallPutFlag, S, X, T, r, b, v)

End Function


'// Forward start options
Public Function ForwardStartOption(CallPutFlag As String, S As Double, alpha As Double, t1 As Double, _
                T As Double, r As Double, b As Double, v As Double) As Double

    ForwardStartOption = S * Exp((b - r) * t1) * GBlackScholes(CallPutFlag, 1, alpha, T - t1, r, b, v)

End Function


'// Time switch options (discrete)
Public Function TimeSwitchOption(CallPutFlag As String, S As Double, X As Double, a As Double, T As Double, m As Integer, dt As Double, r As Double, b As Double, v As Double) As Double
    
    Dim Sum As Double, d As Double
    Dim I As Integer, n As Integer, Z As Integer
    
    n = T / dt
    Sum = 0
    If CallPutFlag = "c" Then
        Z = 1
    ElseIf CallPutFlag = "p" Then
        Z = -1
    End If
    For I = 1 To n
        d = (Log(S / X) + (b - v ^ 2 / 2) * I * dt) / (v * Sqr(I * dt))
        Sum = Sum + CND(Z * d) * dt
    Next
    TimeSwitchOption = a * Exp(-r * T) * Sum + dt * a * Exp(-r * T) * m
End Function


'// Simple chooser options
Public Function SimpleChooser(S As Double, X As Double, t1 As Double, T2 As Double, _
                r As Double, b As Double, v As Double) As Double

    Dim d As Double, y As Double

    d = (Log(S / X) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    y = (Log(S / X) + b * T2 + v ^ 2 * t1 / 2) / (v * Sqr(t1))
  
    SimpleChooser = S * Exp((b - r) * T2) * CND(d) - X * Exp(-r * T2) * CND(d - v * Sqr(T2)) _
    - S * Exp((b - r) * T2) * CND(-y) + X * Exp(-r * T2) * CND(-y + v * Sqr(t1))
End Function


'// Complex chooser options
Public Function ComplexChooser(S As Double, Xc As Double, Xp As Double, T As Double, Tc As Double, _
                Tp As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double, d2 As Double, y1 As Double, y2 As Double
    Dim rho1 As Double, rho2 As Double, I As Double

    I = CriticalValueChooser(S, Xc, Xp, T, Tc, Tp, r, b, v)
    d1 = (Log(S / I) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    y1 = (Log(S / Xc) + (b + v ^ 2 / 2) * Tc) / (v * Sqr(Tc))
    y2 = (Log(S / Xp) + (b + v ^ 2 / 2) * Tp) / (v * Sqr(Tp))
    rho1 = Sqr(T / Tc)
    rho2 = Sqr(T / Tp)
    
    ComplexChooser = S * Exp((b - r) * Tc) * CBND(d1, y1, rho1) - Xc * Exp(-r * Tc) * CBND(d2, y1 - v * Sqr(Tc), rho1) - S * Exp((b - r) * Tp) * CBND(-d1, -y2, rho2) + Xp * Exp(-r * Tp) * CBND(-d2, -y2 + v * Sqr(Tp), rho2)
End Function

'// Critical value complex chooser option
Private Function CriticalValueChooser(S As Double, Xc As Double, Xp As Double, T As Double, _
                Tc As Double, Tp As Double, r As Double, b As Double, v As Double) As Double

    Dim Sv As Double, ci As Double, Pi As Double, epsilon As Double
    Dim dc As Double, dp As Double, yi As Double, di As Double

    Sv = S
    
    ci = GBlackScholes("c", Sv, Xc, Tc - T, r, b, v)
    Pi = GBlackScholes("p", Sv, Xp, Tp - T, r, b, v)
    dc = GDelta("c", Sv, Xc, Tc - T, r, b, v)
    dp = GDelta("p", Sv, Xp, Tp - T, r, b, v)
    yi = ci - Pi
    di = dc - dp
    epsilon = 0.001
    'Newton-Raphson skeprosess
    While Abs(yi) > epsilon
        Sv = Sv - (yi) / di
        ci = GBlackScholes("c", Sv, Xc, Tc - T, r, b, v)
        Pi = GBlackScholes("p", Sv, Xp, Tp - T, r, b, v)
        dc = GDelta("c", Sv, Xc, Tc - T, r, b, v)
        dp = GDelta("p", Sv, Xp, Tp - T, r, b, v)
        yi = ci - Pi
        di = dc - dp
    Wend
    CriticalValueChooser = Sv
End Function


'// Options on options
Public Function OptionsOnOptions(TypeFlag As String, S As Double, X1 As Double, X2 As Double, t1 As Double, _
                T2 As Double, r As Double, b As Double, v As Double) As Double

    Dim y1 As Double, y2 As Double, z1 As Double, z2 As Double
    Dim I As Double, rho As Double, CallPutFlag As String
    
    If TypeFlag = "cc" Or TypeFlag = "pc" Then
        CallPutFlag = "c"
    Else
        CallPutFlag = "p"
    End If
    
    I = CriticalValueOptionsOnOptions(CallPutFlag, X1, X2, T2 - t1, r, b, v)
    
    rho = Sqr(t1 / T2)
    y1 = (Log(S / I) + (b + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    y2 = y1 - v * Sqr(t1)
    z1 = (Log(S / X1) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    z2 = z1 - v * Sqr(T2)

    If TypeFlag = "cc" Then
        OptionsOnOptions = S * Exp((b - r) * T2) * CBND(z1, y1, rho) - X1 * Exp(-r * T2) * CBND(z2, y2, rho) - X2 * Exp(-r * t1) * CND(y2)
    ElseIf TypeFlag = "pc" Then
        OptionsOnOptions = X1 * Exp(-r * T2) * CBND(z2, -y2, -rho) - S * Exp((b - r) * T2) * CBND(z1, -y1, -rho) + X2 * Exp(-r * t1) * CND(-y2)
    ElseIf TypeFlag = "cp" Then
        OptionsOnOptions = X1 * Exp(-r * T2) * CBND(-z2, -y2, rho) - S * Exp((b - r) * T2) * CBND(-z1, -y1, rho) - X2 * Exp(-r * t1) * CND(-y2)
    ElseIf TypeFlag = "pp" Then
        OptionsOnOptions = S * Exp((b - r) * T2) * CBND(-z1, y1, -rho) - X1 * Exp(-r * T2) * CBND(-z2, y2, -rho) + Exp(-r * t1) * X2 * CND(y2)
    End If
End Function
'// Calculation of critical price options on options
Private Function CriticalValueOptionsOnOptions(CallPutFlag As String, X1 As Double, X2 As Double, T As Double, _
                r As Double, b As Double, v As Double) As Double

    Dim Si As Double, ci As Double, di As Double, epsilon As Double
    
    Si = X1
    ci = GBlackScholes(CallPutFlag, Si, X1, T, r, b, v)
    di = GDelta(CallPutFlag, Si, X1, T, r, b, v)
    epsilon = 0.000001
    '// Newton-Raphson algorithm
    While Abs(ci - X2) > epsilon
        Si = Si - (ci - X2) / di
        ci = GBlackScholes(CallPutFlag, Si, X1, T, r, b, v)
        di = GDelta(CallPutFlag, Si, X1, T, r, b, v)
    Wend
    CriticalValueOptionsOnOptions = Si
End Function


'// Writer extendible options
Public Function ExtendibleWriter(CallPutFlag As String, S As Double, X1 As Double, X2 As Double, t1 As Double, _
                T2 As Double, r As Double, b As Double, v As Double) As Double

    Dim rho As Double, z1 As Double, z2 As Double
    rho = Sqr(t1 / T2)
    z1 = (Log(S / X2) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    z2 = (Log(S / X1) + (b + v ^ 2 / 2) * t1) / (v * Sqr(t1))

    If CallPutFlag = "c" Then
        ExtendibleWriter = GBlackScholes(CallPutFlag, S, X1, t1, r, b, v) + S * Exp((b - r) * T2) * CBND(z1, -z2, -rho) - X2 * Exp(-r * T2) * CBND(z1 - Sqr(v ^ 2 * T2), -z2 + Sqr(v ^ 2 * t1), -rho)
    ElseIf CallPutFlag = "p" Then
        ExtendibleWriter = GBlackScholes(CallPutFlag, S, X1, t1, r, b, v) + X2 * Exp(-r * T2) * CBND(-z1 + Sqr(v ^ 2 * T2), z2 - Sqr(v ^ 2 * t1), -rho) - S * Exp((b - r) * T2) * CBND(-z1, z2, -rho)
    End If
End Function


'// Two asset correlation options
Public Function TwoAssetCorrelation(CallPutFlag As String, S1 As Double, S2 As Double, X1 As Double, X2 As Double, T As Double, _
                b1 As Double, b2 As Double, r As Double, v1 As Double, v2 As Double, rho As Double)

    Dim y1 As Double, y2 As Double
   
    y1 = (Log(S1 / X1) + (b1 - v1 ^ 2 / 2) * T) / (v1 * Sqr(T))
    y2 = (Log(S2 / X2) + (b2 - v2 ^ 2 / 2) * T) / (v2 * Sqr(T))
    
    If CallPutFlag = "c" Then
        TwoAssetCorrelation = S2 * Exp((b2 - r) * T) * CBND(y2 + v2 * Sqr(T), y1 + rho * v2 * Sqr(T), rho) _
        - X2 * Exp(-r * T) * CBND(y2, y1, rho)
    ElseIf CallPutFlag = "p" Then
         TwoAssetCorrelation = X2 * Exp(-r * T) * CBND(-y2, -y1, rho) _
         - S2 * Exp((b2 - r) * T) * CBND(-y2 - v2 * Sqr(T), -y1 - rho * v2 * Sqr(T), rho)
    End If
End Function


'// European option to exchange one asset for another
Public Function EuropeanExchangeOption(S1 As Double, S2 As Double, Q1 As Double, Q2 As Double, T As Double, r As Double, b1 As Double, _
                b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double

    Dim v As Double, d1 As Double, d2 As Double

    v = Sqr(v1 ^ 2 + v2 ^ 2 - 2 * rho * v1 * v2)
    d1 = (Log(Q1 * S1 / (Q2 * S2)) + (b1 - b2 + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
 
    EuropeanExchangeOption = Q1 * S1 * Exp((b1 - r) * T) * CND(d1) - Q2 * S2 * Exp((b2 - r) * T) * CND(d2)
End Function


'// American option to exchange one asset for another
Public Function AmericanExchangeOption(S1 As Double, S2 As Double, Q1 As Double, Q2 As Double, T As Double, _
            r As Double, b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double
    Dim v As Double
    v = Sqr(v1 ^ 2 + v2 ^ 2 - 2 * rho * v1 * v2)
    AmericanExchangeOption = BSAmericanApprox("c", Q1 * S1, Q2 * S2, T, r - b2, b1 - b2, v)
End Function


'// Exchange options on exchange options
Public Function ExchangeExchangeOption(TypeFlag As Integer, S1 As Double, S2 As Double, q As Double, t1 As Double, T2 As Double, r As Double, b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double
    
    Dim I As Double, I1 As Double
    Dim d1 As Double, d2 As Double
    Dim d3 As Double, d4 As Double
    Dim y1 As Double, y2 As Double
    Dim y3 As Double, y4 As Double
    Dim v As Double, id As Integer
    
    v = Sqr(v1 ^ 2 + v2 ^ 2 - 2 * rho * v1 * v2)
    I1 = S1 * Exp((b1 - r) * (T2 - t1)) / (S2 * Exp((b2 - r) * (T2 - t1)))
    
    If TypeFlag = 1 Or TypeFlag = 2 Then
        id = 1
    Else
        id = 2
    End If
    
    I = CriticalPrice(id, I1, t1, T2, v, q)
    d1 = (Log(S1 / (I * S2)) + (b1 - b2 + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    d2 = d1 - v * Sqr(t1)
    d3 = (Log((I * S2) / S1) + (b2 - b1 + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    d4 = d3 - v * Sqr(t1)
    y1 = (Log(S1 / S2) + (b1 - b2 + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    y2 = y1 - v * Sqr(T2)
    y3 = (Log(S2 / S1) + (b2 - b1 + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    y4 = y3 - v * Sqr(T2)
    
    If TypeFlag = 1 Then
        ExchangeExchangeOption = -S2 * Exp((b2 - r) * T2) * CBND(d2, y2, Sqr(t1 / T2)) + S1 * Exp((b1 - r) * T2) * CBND(d1, y1, Sqr(t1 / T2)) - q * S2 * Exp((b2 - r) * t1) * CND(d2)
    ElseIf TypeFlag = 2 Then
        ExchangeExchangeOption = S2 * Exp((b2 - r) * T2) * CBND(d3, y2, -Sqr(t1 / T2)) - S1 * Exp((b1 - r) * T2) * CBND(d4, y1, -Sqr(t1 / T2)) + q * S2 * Exp((b2 - r) * t1) * CND(d3)
    ElseIf TypeFlag = 3 Then
        ExchangeExchangeOption = S2 * Exp((b2 - r) * T2) * CBND(d3, y3, Sqr(t1 / T2)) - S1 * Exp((b1 - r) * T2) * CBND(d4, y4, Sqr(t1 / T2)) - q * S2 * Exp((b2 - r) * t1) * CND(d3)
    ElseIf TypeFlag = 4 Then
        ExchangeExchangeOption = -S2 * Exp((b2 - r) * T2) * CBND(d2, y3, -Sqr(t1 / T2)) + S1 * Exp((b1 - r) * T2) * CBND(d1, y4, -Sqr(t1 / T2)) + q * S2 * Exp((b2 - r) * t1) * CND(d2)
    End If
End Function

'// Numerical search algorithm to find critical price I
Private Function CriticalPrice(id As Integer, I1 As Double, t1 As Double, T2 As Double, v As Double, q As Double) As Double
    Dim Ii As Double, yi As Double, di As Double
    Dim epsilon As Double
    
    Ii = I1
    yi = CriticalPart3(id, Ii, t1, T2, v)
    di = CriticalPart2(id, Ii, t1, T2, v)
    epsilon = 0.00001
    While Abs(yi - q) > epsilon
        Ii = Ii - (yi - q) / di
        yi = CriticalPart3(id, Ii, t1, T2, v)
        di = CriticalPart2(id, Ii, t1, T2, v)
    Wend
    CriticalPrice = Ii
End Function
Private Function CriticalPart2(id As Integer, I As Double, t1 As Double, T2 As Double, v As Double) As Double
    Dim z1 As Double, z2 As Double
    If id = 1 Then
        z1 = (Log(I) + v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        CriticalPart2 = CND(z1)
    ElseIf id = 2 Then
        z2 = (-Log(I) - v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        CriticalPart2 = -CND(z2)
    End If
End Function
Private Function CriticalPart3(id As Integer, I As Double, t1 As Double, T2 As Double, v As Double) As Double
    Dim z1 As Double, z2 As Double
    If id = 1 Then
        z1 = (Log(I) + v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        z2 = (Log(I) - v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        CriticalPart3 = I * CND(z1) - CND(z2)
    ElseIf id = 2 Then
        z1 = (-Log(I) + v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        z2 = (-Log(I) - v ^ 2 / 2 * (T2 - t1)) / (v * Sqr(T2 - t1))
        CriticalPart3 = CND(z1) - I * CND(z2)
    End If
End Function


'// Options on the maximum or minimum of two risky assets
Public Function OptionsOnTheMaxMin(TypeFlag As String, S1 As Double, S2 As Double, X As Double, T As Double, r As Double, _
        b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double

    Dim v As Double, rho1 As Double, rho2 As Double
    Dim d As Double, y1 As Double, y2 As Double
    
    v = Sqr(v1 ^ 2 + v2 ^ 2 - 2 * rho * v1 * v2)
    rho1 = (v1 - rho * v2) / v
    rho2 = (v2 - rho * v1) / v
    d = (Log(S1 / S2) + (b1 - b2 + v ^ 2 / 2) * T) / (v * Sqr(T))
    y1 = (Log(S1 / X) + (b1 + v1 ^ 2 / 2) * T) / (v1 * Sqr(T))
    y2 = (Log(S2 / X) + (b2 + v2 ^ 2 / 2) * T) / (v2 * Sqr(T))
  
    If TypeFlag = "cmin" Then
        OptionsOnTheMaxMin = S1 * Exp((b1 - r) * T) * CBND(y1, -d, -rho1) + S2 * Exp((b2 - r) * T) * CBND(y2, d - v * Sqr(T), -rho2) - X * Exp(-r * T) * CBND(y1 - v1 * Sqr(T), y2 - v2 * Sqr(T), rho)
    ElseIf TypeFlag = "cmax" Then
        OptionsOnTheMaxMin = S1 * Exp((b1 - r) * T) * CBND(y1, d, rho1) + S2 * Exp((b2 - r) * T) * CBND(y2, -d + v * Sqr(T), rho2) - X * Exp(-r * T) * (1 - CBND(-y1 + v1 * Sqr(T), -y2 + v2 * Sqr(T), rho))
    ElseIf TypeFlag = "pmin" Then
        OptionsOnTheMaxMin = X * Exp(-r * T) - S1 * Exp((b1 - r) * T) + EuropeanExchangeOption(S1, S2, 1, 1, T, r, b1, b2, v1, v2, rho) + OptionsOnTheMaxMin("cmin", S1, S2, X, T, r, b1, b2, v1, v2, rho)
    ElseIf TypeFlag = "pmax" Then
        OptionsOnTheMaxMin = X * Exp(-r * T) - S2 * Exp((b2 - r) * T) - EuropeanExchangeOption(S1, S2, 1, 1, T, r, b1, b2, v1, v2, rho) + OptionsOnTheMaxMin("cmax", S1, S2, X, T, r, b1, b2, v1, v2, rho)
    End If
End Function


'// Spread option approximation
Function SpreadApproximation(CallPutFlag As String, f1 As Double, f2 As Double, X As Double, T As Double, _
                r As Double, v1 As Double, v2 As Double, rho As Double) As Double

    Dim v As Double, F As Double
    Dim d1 As Double, d2 As Double
   
    v = Sqr(v1 ^ 2 + (v2 * f2 / (f2 + X)) ^ 2 - 2 * rho * v1 * v2 * f2 / (f2 + X))
    F = f1 / (f2 + X)
    
    SpreadApproximation = GBlackScholes(CallPutFlag, F, 1, T, r, 0, v) * (f2 + X)
End Function


'// Floating strike lookback options
Function FloatingStrikeLookback(CallPutFlag As String, S As Double, SMin As Double, SMax As Double, T As Double, _
            r As Double, b As Double, v As Double) As Double

    Dim a1 As Double, a2 As Double, m As Double
    
    If CallPutFlag = "c" Then
        m = SMin
    ElseIf CallPutFlag = "p" Then
        m = SMax
    End If
    
     a1 = (Log(S / m) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
     a2 = a1 - v * Sqr(T)

    If CallPutFlag = "c" Then
        FloatingStrikeLookback = S * Exp((b - r) * T) * CND(a1) - m * Exp(-r * T) * CND(a2) + _
        Exp(-r * T) * v ^ 2 / (2 * b) * S * ((S / m) ^ (-2 * b / v ^ 2) * CND(-a1 + 2 * b / v * Sqr(T)) - Exp(b * T) * CND(-a1))
    ElseIf CallPutFlag = "p" Then
        FloatingStrikeLookback = m * Exp(-r * T) * CND(-a2) - S * Exp((b - r) * T) * CND(-a1) + _
        Exp(-r * T) * v ^ 2 / (2 * b) * S * (-(S / m) ^ (-2 * b / v ^ 2) * CND(a1 - 2 * b / v * Sqr(T)) + Exp(b * T) * CND(a1))
    End If
End Function


'// Fixed strike lookback options
Public Function FixedStrikeLookback(CallPutFlag As String, S As Double, SMin As Double, SMax As Double, X As Double, _
                 T As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double, d2 As Double
    Dim e1 As Double, e2 As Double, m As Double
    
    If CallPutFlag = "c" Then
        m = SMax
    ElseIf CallPutFlag = "p" Then
        m = SMin
    End If
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    e1 = (Log(S / m) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    e2 = e1 - v * Sqr(T)
    
    If CallPutFlag = "c" And X > m Then
        FixedStrikeLookback = S * Exp((b - r) * T) * CND(d1) - X * Exp(-r * T) * CND(d2) _
        + S * Exp(-r * T) * v ^ 2 / (2 * b) * (-(S / X) ^ (-2 * b / v ^ 2) * CND(d1 - 2 * b / v * Sqr(T)) + Exp(b * T) * CND(d1))
    ElseIf CallPutFlag = "c" And X <= m Then
        FixedStrikeLookback = Exp(-r * T) * (m - X) + S * Exp((b - r) * T) * CND(e1) - Exp(-r * T) * m * CND(e2) _
        + S * Exp(-r * T) * v ^ 2 / (2 * b) * (-(S / m) ^ (-2 * b / v ^ 2) * CND(e1 - 2 * b / v * Sqr(T)) + Exp(b * T) * CND(e1))
    ElseIf CallPutFlag = "p" And X < m Then
        FixedStrikeLookback = -S * Exp((b - r) * T) * CND(-d1) + X * Exp(-r * T) * CND(-d1 + v * Sqr(T)) _
        + S * Exp(-r * T) * v ^ 2 / (2 * b) * ((S / X) ^ (-2 * b / v ^ 2) * CND(-d1 + 2 * b / v * Sqr(T)) - Exp(b * T) * CND(-d1))
    ElseIf CallPutFlag = "p" And X >= m Then
        FixedStrikeLookback = Exp(-r * T) * (X - m) - S * Exp((b - r) * T) * CND(-e1) + Exp(-r * T) * m * CND(-e1 + v * Sqr(T)) _
        + Exp(-r * T) * v ^ 2 / (2 * b) * S * ((S / m) ^ (-2 * b / v ^ 2) * CND(-e1 + 2 * b / v * Sqr(T)) - Exp(b * T) * CND(-e1))
    End If
End Function


'// Partial-time floating strike lookback options
Public Function PartialFloatLB(CallPutFlag As String, S As Double, SMin As Double, SMax As Double, t1 As Double, _
                T2 As Double, r As Double, b As Double, v As Double, lambda As Double)
   
    Dim d1 As Double, d2 As Double
    Dim e1 As Double, e2 As Double
    Dim f1 As Double, f2 As Double
    Dim g1 As Double, g2 As Double, m As Double
    Dim part1 As Double, part2 As Double, part3 As Double
    
    If CallPutFlag = "c" Then
        m = SMin
    ElseIf CallPutFlag = "p" Then
        m = SMax
    End If
    
    d1 = (Log(S / m) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    d2 = d1 - v * Sqr(T2)
    e1 = (b + v ^ 2 / 2) * (T2 - t1) / (v * Sqr(T2 - t1))
    e2 = e1 - v * Sqr(T2 - t1)
    f1 = (Log(S / m) + (b + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    f2 = f1 - v * Sqr(t1)
    g1 = Log(lambda) / (v * Sqr(T2))
    g2 = Log(lambda) / (v * Sqr(T2 - t1))

    If CallPutFlag = "c" Then
        part1 = S * Exp((b - r) * T2) * CND(d1 - g1) - lambda * m * Exp(-r * T2) * CND(d2 - g1)
        part2 = Exp(-r * T2) * v ^ 2 / (2 * b) * lambda * S * ((S / m) ^ (-2 * b / v ^ 2) * CBND(-f1 + 2 * b * Sqr(t1) / v, -d1 + 2 * b * Sqr(T2) / v - g1, Sqr(t1 / T2)) _
        - Exp(b * T2) * lambda ^ (2 * b / v ^ 2) * CBND(-d1 - g1, e1 + g2, -Sqr(1 - t1 / T2))) _
        + S * Exp((b - r) * T2) * CBND(-d1 + g1, e1 - g2, -Sqr(1 - t1 / T2))
        part3 = Exp(-r * T2) * lambda * m * CBND(-f2, d2 - g1, -Sqr(t1 / T2)) _
        - Exp(-b * (T2 - t1)) * Exp((b - r) * T2) * (1 + v ^ 2 / (2 * b)) * lambda * S * CND(e2 - g2) * CND(-f1)
    
    ElseIf CallPutFlag = "p" Then
        part1 = lambda * m * Exp(-r * T2) * CND(-d2 + g1) - S * Exp((b - r) * T2) * CND(-d1 + g1)
        part2 = -Exp(-r * T2) * v ^ 2 / (2 * b) * lambda * S * ((S / m) ^ (-2 * b / v ^ 2) * CBND(f1 - 2 * b * Sqr(t1) / v, d1 - 2 * b * Sqr(T2) / v + g1, Sqr(t1 / T2)) _
        - Exp(b * T2) * lambda ^ (2 * b / v ^ 2) * CBND(d1 + g1, -e1 - g2, -Sqr(1 - t1 / T2))) _
        - S * Exp((b - r) * T2) * CBND(d1 - g1, -e1 + g2, -Sqr(1 - t1 / T2))
        part3 = -Exp(-r * T2) * lambda * m * CBND(f2, -d2 + g1, -Sqr(t1 / T2)) _
        + Exp(-b * (T2 - t1)) * Exp((b - r) * T2) * (1 + v ^ 2 / (2 * b)) * lambda * S * CND(-e2 + g2) * CND(f1)
  End If
  PartialFloatLB = part1 + part2 + part3
End Function


'// Partial-time fixed strike lookback options
Public Function PartialFixedLB(CallPutFlag As String, S As Double, X As Double, t1 As Double, _
                T2 As Double, r As Double, b As Double, v As Double) As Double

    Dim d1 As Double, d2 As Double
    Dim e1 As Double, e2 As Double
    Dim f1 As Double, f2 As Double

    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    d2 = d1 - v * Sqr(T2)
    e1 = ((b + v ^ 2 / 2) * (T2 - t1)) / (v * Sqr(T2 - t1))
    e2 = e1 - v * Sqr(T2 - t1)
    f1 = (Log(S / X) + (b + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    f2 = f1 - v * Sqr(t1)
    If CallPutFlag = "c" Then
        PartialFixedLB = S * Exp((b - r) * T2) * CND(d1) - Exp(-r * T2) * X * CND(d2) + S * Exp(-r * T2) * v ^ 2 / (2 * b) * (-(S / X) ^ (-2 * b / v ^ 2) * CBND(d1 - 2 * b * Sqr(T2) / v, -f1 + 2 * b * Sqr(t1) / v, -Sqr(t1 / T2)) + Exp(b * T2) * CBND(e1, d1, Sqr(1 - t1 / T2))) - S * Exp((b - r) * T2) * CBND(-e1, d1, -Sqr(1 - t1 / T2)) - X * Exp(-r * T2) * CBND(f2, -d2, -Sqr(t1 / T2)) + Exp(-b * (T2 - t1)) * (1 - v ^ 2 / (2 * b)) * S * Exp((b - r) * T2) * CND(f1) * CND(-e2)
    ElseIf CallPutFlag = "p" Then
        PartialFixedLB = X * Exp(-r * T2) * CND(-d2) - S * Exp((b - r) * T2) * CND(-d1) + S * Exp(-r * T2) * v ^ 2 / (2 * b) * ((S / X) ^ (-2 * b / v ^ 2) * CBND(-d1 + 2 * b * Sqr(T2) / v, f1 - 2 * b * Sqr(t1) / v, -Sqr(t1 / T2)) - Exp(b * T2) * CBND(-e1, -d1, Sqr(1 - t1 / T2))) + S * Exp((b - r) * T2) * CBND(e1, -d1, -Sqr(1 - t1 / T2)) + X * Exp(-r * T2) * CBND(-f2, d2, -Sqr(t1 / T2)) - Exp(-b * (T2 - t1)) * (1 - v ^ 2 / (2 * b)) * S * Exp((b - r) * T2) * CND(-f1) * CND(e2)
    End If
End Function


'// Extreme spread options
Public Function ExtremeSpreadOption(TypeFlag As Integer, S As Double, SMin As Double, SMax As Double, t1 As Double, T As Double, r As Double, b As Double, v As Double) As Double
            
        Dim m As Double, Mo As Double
        Dim mu1 As Double, mu As Double
        Dim kappa As Integer, eta As Integer
        
        If TypeFlag = 1 Or TypeFlag = 3 Then
            eta = 1
        Else
            eta = -1
        End If
        If TypeFlag = 1 Or TypeFlag = 2 Then
            kappa = 1
        Else
            kappa = -1
        End If
            
        If kappa * eta = 1 Then
            Mo = SMax
        ElseIf kappa * eta = -1 Then
            Mo = SMin
        End If
        
        mu1 = b - v ^ 2 / 2
        mu = mu1 + v ^ 2
        m = Log(Mo / S)
        If kappa = 1 Then '// Extreme Spread Option
            ExtremeSpreadOption = eta * (S * Exp((b - r) * T) * (1 + v ^ 2 / (2 * b)) * CND(eta * (-m + mu * T) / (v * Sqr(T))) - Exp(-r * (T - t1)) * S * Exp((b - r) * T) * (1 + v ^ 2 / (2 * b)) * CND(eta * (-m + mu * t1) / (v * Sqr(t1))) _
            + Exp(-r * T) * Mo * CND(eta * (m - mu1 * T) / (v * Sqr(T))) - Exp(-r * T) * Mo * v ^ 2 / (2 * b) * Exp(2 * mu1 * m / v ^ 2) * CND(eta * (-m - mu1 * T) / (v * Sqr(T))) _
            - Exp(-r * T) * Mo * CND(eta * (m - mu1 * t1) / (v * Sqr(t1))) + Exp(-r * T) * Mo * v ^ 2 / (2 * b) * Exp(2 * mu1 * m / v ^ 2) * CND(eta * (-m - mu1 * t1) / (v * Sqr(t1))))
        ElseIf kappa = -1 Then  '// Reverse Extreme Spread Option
            ExtremeSpreadOption = -eta * (S * Exp((b - r) * T) * (1 + v ^ 2 / (2 * b)) * CND(eta * (m - mu * T) / (v * Sqr(T))) + Exp(-r * T) * Mo * CND(eta * (-m + mu1 * T) / (v * Sqr(T))) _
            - Exp(-r * T) * Mo * v ^ 2 / (2 * b) * Exp(2 * mu1 * m / v ^ 2) * CND(eta * (m + mu1 * T) / (v * Sqr(T))) - S * Exp((b - r) * T) * (1 + v ^ 2 / (2 * b)) * CND(eta * (-mu * (T - t1)) / (v * Sqr(T - t1))) _
            - Exp(-r * (T - t1)) * S * Exp((b - r) * T) * (1 - v ^ 2 / (2 * b)) * CND(eta * (mu1 * (T - t1)) / (v * Sqr(T - t1))))
        End If
End Function


'// Standard barrier options
Function StandardBarrier(TypeFlag As String, S As Double, X As Double, H As Double, K As Double, T As Double, _
            r As Double, b As Double, v As Double)
    
    Dim mu As Double
    Dim lambda As Double
    Dim X1 As Double, X2 As Double
    Dim y1 As Double, y2 As Double
    Dim Z As Double
    
    Dim eta As Integer    'Binary variable that can take the value of 1 or -1
    Dim phi As Integer    'Binary variable that can take the value of 1 or -1
    
    Dim f1 As Double    'Equal to formula "A" in the book
    Dim f2 As Double    'Equal to formula "B" in the book
    Dim f3 As Double    'Equal to formula "C" in the book
    Dim f4 As Double    'Equal to formula "D" in the book
    Dim f5 As Double    'Equal to formula "E" in the book
    Dim f6 As Double    'Equal to formula "F" in the book

    mu = (b - v ^ 2 / 2) / v ^ 2
    lambda = Sqr(mu ^ 2 + 2 * r / v ^ 2)
    X1 = Log(S / X) / (v * Sqr(T)) + (1 + mu) * v * Sqr(T)
    X2 = Log(S / H) / (v * Sqr(T)) + (1 + mu) * v * Sqr(T)
    y1 = Log(H ^ 2 / (S * X)) / (v * Sqr(T)) + (1 + mu) * v * Sqr(T)
    y2 = Log(H / S) / (v * Sqr(T)) + (1 + mu) * v * Sqr(T)
    Z = Log(H / S) / (v * Sqr(T)) + lambda * v * Sqr(T)
    
    If TypeFlag = "cdi" Or TypeFlag = "cdo" Then
        eta = 1
        phi = 1
    ElseIf TypeFlag = "cui" Or TypeFlag = "cuo" Then
        eta = -1
        phi = 1
    ElseIf TypeFlag = "pdi" Or TypeFlag = "pdo" Then
        eta = 1
        phi = -1
    ElseIf TypeFlag = "pui" Or TypeFlag = "puo" Then
        eta = -1
        phi = -1
    End If
    
    f1 = phi * S * Exp((b - r) * T) * CND(phi * X1) - phi * X * Exp(-r * T) * CND(phi * X1 - phi * v * Sqr(T))
    f2 = phi * S * Exp((b - r) * T) * CND(phi * X2) - phi * X * Exp(-r * T) * CND(phi * X2 - phi * v * Sqr(T))
    f3 = phi * S * Exp((b - r) * T) * (H / S) ^ (2 * (mu + 1)) * CND(eta * y1) - phi * X * Exp(-r * T) * (H / S) ^ (2 * mu) * CND(eta * y1 - eta * v * Sqr(T))
    f4 = phi * S * Exp((b - r) * T) * (H / S) ^ (2 * (mu + 1)) * CND(eta * y2) - phi * X * Exp(-r * T) * (H / S) ^ (2 * mu) * CND(eta * y2 - eta * v * Sqr(T))
    f5 = K * Exp(-r * T) * (CND(eta * X2 - eta * v * Sqr(T)) - (H / S) ^ (2 * mu) * CND(eta * y2 - eta * v * Sqr(T)))
    f6 = K * ((H / S) ^ (mu + lambda) * CND(eta * Z) + (H / S) ^ (mu - lambda) * CND(eta * Z - 2 * eta * lambda * v * Sqr(T)))
    
    
    If X > H Then
        Select Case TypeFlag
            Case Is = "cdi"
                StandardBarrier = f3 + f5
            Case Is = "cui"
                StandardBarrier = f1 + f5
            Case Is = "pdi"
                StandardBarrier = f2 - f3 + f4 + f5
            Case Is = "pui"
                StandardBarrier = f1 - f2 + f4 + f5
            Case Is = "cdo"
                StandardBarrier = f1 - f3 + f6
            Case Is = "cuo"
                StandardBarrier = f6
            Case Is = "pdo"
                StandardBarrier = f1 - f2 + f3 - f4 + f6
            Case Is = "puo"
                StandardBarrier = f2 - f4 + f6
            End Select
    ElseIf X < H Then
        Select Case TypeFlag
            Case Is = "cdi"
                StandardBarrier = f1 - f2 + f4 + f5
            Case Is = "cui"
                StandardBarrier = f2 - f3 + f4 + f5
            Case Is = "pdi"
                StandardBarrier = f1 + f5
            Case Is = "pui"
                StandardBarrier = f3 + f5
            Case Is = "cdo"
                StandardBarrier = f2 + f6 - f4
            Case Is = "cuo"
                StandardBarrier = f1 - f2 + f3 - f4 + f6
            Case Is = "pdo"
                StandardBarrier = f6
            Case Is = "puo"
                StandardBarrier = f1 - f3 + f6
        End Select
    End If
End Function


'// Double barrier options
Function DoubleBarrier(TypeFlag As String, S As Double, X As Double, L As Double, U As Double, T As Double, _
        r As Double, b As Double, v As Double, delta1 As Double, delta2 As Double) As Double
    
    Dim E As Double, F As Double
    Dim Sum1 As Double, Sum2 As Double
    Dim d1 As Double, d2 As Double
    Dim d3 As Double, d4 As Double
    Dim mu1 As Double, mu2 As Double, mu3 As Double
    Dim OutValue As Double, n As Integer
    
    F = U * Exp(delta1 * T)
    E = L * Exp(delta1 * T)
    Sum1 = 0
    Sum2 = 0
    If TypeFlag = "co" Or TypeFlag = "ci" Then
        For n = -5 To 5
            d1 = (Log(S * U ^ (2 * n) / (X * L ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d2 = (Log(S * U ^ (2 * n) / (F * L ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d3 = (Log(L ^ (2 * n + 2) / (X * S * U ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d4 = (Log(L ^ (2 * n + 2) / (F * S * U ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ^ 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ^ 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ^ 2 + 1
            Sum1 = Sum1 + (U ^ n / L ^ n) ^ mu1 * (L / S) ^ mu2 * (CND(d1) - CND(d2)) - (L ^ (n + 1) / (U ^ n * S)) ^ mu3 * (CND(d3) - CND(d4))
            Sum2 = Sum2 + (U ^ n / L ^ n) ^ (mu1 - 2) * (L / S) ^ mu2 * (CND(d1 - v * Sqr(T)) - CND(d2 - v * Sqr(T))) - (L ^ (n + 1) / (U ^ n * S)) ^ (mu3 - 2) * (CND(d3 - v * Sqr(T)) - CND(d4 - v * Sqr(T)))
        Next
        OutValue = S * Exp((b - r) * T) * Sum1 - X * Exp(-r * T) * Sum2
    ElseIf TypeFlag = "po" Or TypeFlag = "pi" Then
        For n = -5 To 5
            d1 = (Log(S * U ^ (2 * n) / (E * L ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d2 = (Log(S * U ^ (2 * n) / (X * L ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d3 = (Log(L ^ (2 * n + 2) / (E * S * U ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            d4 = (Log(L ^ (2 * n + 2) / (X * S * U ^ (2 * n))) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ^ 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ^ 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ^ 2 + 1
            Sum1 = Sum1 + (U ^ n / L ^ n) ^ mu1 * (L / S) ^ mu2 * (CND(d1) - CND(d2)) - (L ^ (n + 1) / (U ^ n * S)) ^ mu3 * (CND(d3) - CND(d4))
            Sum2 = Sum2 + (U ^ n / L ^ n) ^ (mu1 - 2) * (L / S) ^ mu2 * (CND(d1 - v * Sqr(T)) - CND(d2 - v * Sqr(T))) - (L ^ (n + 1) / (U ^ n * S)) ^ (mu3 - 2) * (CND(d3 - v * Sqr(T)) - CND(d4 - v * Sqr(T)))
        Next
        OutValue = X * Exp(-r * T) * Sum2 - S * Exp((b - r) * T) * Sum1
    End If
    If TypeFlag = "co" Or TypeFlag = "po" Then
        DoubleBarrier = OutValue
    ElseIf TypeFlag = "ci" Then
        DoubleBarrier = GBlackScholes("c", S, X, T, r, b, v) - OutValue
    ElseIf TypeFlag = "pi" Then
        DoubleBarrier = GBlackScholes("p", S, X, T, r, b, v) - OutValue
    End If
End Function


'// Partial-time singel asset barrier options
Public Function PartialTimeBarrier(TypeFlag As String, S As Double, X As Double, H As Double, _
                t1 As Double, T2 As Double, r As Double, b As Double, v As Double) As Double
    
    Dim d1 As Double, d2 As Double
    Dim f1 As Double, f2 As Double
    Dim e1 As Double, e2 As Double
    Dim e3 As Double, e4 As Double
    Dim g1 As Double, g2 As Double
    Dim g3 As Double, g4 As Double
    Dim mu As Double, rho As Double, eta As Integer
    Dim z1 As Double, z2 As Double, z3 As Double
    Dim z4 As Double, z5 As Double, z6 As Double
    Dim z7 As Double, z8 As Double
    
    If TypeFlag = "cdoA" Then
        eta = 1
    ElseIf TypeFlag = "cuoA" Then
        eta = -1
    End If
    
    d1 = (Log(S / X) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    d2 = d1 - v * Sqr(T2)
    f1 = (Log(S / X) + 2 * Log(H / S) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    f2 = f1 - v * Sqr(T2)
    e1 = (Log(S / H) + (b + v ^ 2 / 2) * t1) / (v * Sqr(t1))
    e2 = e1 - v * Sqr(t1)
    e3 = e1 + 2 * Log(H / S) / (v * Sqr(t1))
    e4 = e3 - v * Sqr(t1)
    mu = (b - v ^ 2 / 2) / v ^ 2
    rho = Sqr(t1 / T2)
    g1 = (Log(S / H) + (b + v ^ 2 / 2) * T2) / (v * Sqr(T2))
    g2 = g1 - v * Sqr(T2)
    g3 = g1 + 2 * Log(H / S) / (v * Sqr(T2))
    g4 = g3 - v * Sqr(T2)
    
    z1 = CND(e2) - (H / S) ^ (2 * mu) * CND(e4)
    z2 = CND(-e2) - (H / S) ^ (2 * mu) * CND(-e4)
    z3 = CBND(g2, e2, rho) - (H / S) ^ (2 * mu) * CBND(g4, -e4, -rho)
    z4 = CBND(-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4, e4, -rho)
    z5 = CND(e1) - (H / S) ^ (2 * (mu + 1)) * CND(e3)
    z6 = CND(-e1) - (H / S) ^ (2 * (mu + 1)) * CND(-e3)
    z7 = CBND(g1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(g3, -e3, -rho)
    z8 = CBND(-g1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(-g3, e3, -rho)
    
    If TypeFlag = "cdoA" Or TypeFlag = "cuoA" Then '// call down-and out and up-and-out type A
        PartialTimeBarrier = S * Exp((b - r) * T2) * (CBND(d1, eta * e1, eta * rho) - (H / S) ^ (2 * (mu + 1)) * CBND(f1, eta * e3, eta * rho)) _
        - X * Exp(-r * T2) * (CBND(d2, eta * e2, eta * rho) - (H / S) ^ (2 * mu) * CBND(f2, eta * e4, eta * rho))
    ElseIf TypeFlag = "cdoB2" And X < H Then  '// call down-and-out type B2
        PartialTimeBarrier = S * Exp((b - r) * T2) * (CBND(g1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(g3, -e3, -rho)) _
        - X * Exp(-r * T2) * (CBND(g2, e2, rho) - (H / S) ^ (2 * mu) * CBND(g4, -e4, -rho))
    ElseIf TypeFlag = "cdoB2" And X > H Then
        PartialTimeBarrier = PartialTimeBarrier("coB1", S, X, H, t1, T2, r, b, v)
    ElseIf TypeFlag = "cuoB2" And X < H Then  '// call up-and-out type B2
        PartialTimeBarrier = S * Exp((b - r) * T2) * (CBND(-g1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(-g3, e3, -rho)) _
        - X * Exp(-r * T2) * (CBND(-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4, e4, -rho)) _
        - S * Exp((b - r) * T2) * (CBND(-d1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(e3, -f1, -rho)) _
        + X * Exp(-r * T2) * (CBND(-d2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(e4, -f2, -rho))
    ElseIf TypeFlag = "coB1" And X > H Then  '// call out type B1
        PartialTimeBarrier = S * Exp((b - r) * T2) * (CBND(d1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(f1, -e3, -rho)) _
        - X * Exp(-r * T2) * (CBND(d2, e2, rho) - (H / S) ^ (2 * mu) * CBND(f2, -e4, -rho))
    ElseIf TypeFlag = "coB1" And X < H Then
        PartialTimeBarrier = S * Exp((b - r) * T2) * (CBND(-g1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(-g3, e3, -rho)) _
        - X * Exp(-r * T2) * (CBND(-g2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-g4, e4, -rho)) _
        - S * Exp((b - r) * T2) * (CBND(-d1, -e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(-f1, e3, -rho)) _
        + X * Exp(-r * T2) * (CBND(-d2, -e2, rho) - (H / S) ^ (2 * mu) * CBND(-f2, e4, -rho)) _
        + S * Exp((b - r) * T2) * (CBND(g1, e1, rho) - (H / S) ^ (2 * (mu + 1)) * CBND(g3, -e3, -rho)) _
        - X * Exp(-r * T2) * (CBND(g2, e2, rho) - (H / S) ^ (2 * mu) * CBND(g4, -e4, -rho))
    ElseIf TypeFlag = "pdoA" Then  '// put down-and out and up-and-out type A
        PartialTimeBarrier = PartialTimeBarrier("cdoA", S, X, H, t1, T2, r, b, v) - S * Exp((b - r) * T2) * z5 + X * Exp(-r * T2) * z1
    ElseIf TypeFlag = "puoA" Then
        PartialTimeBarrier = PartialTimeBarrier("cuoA", S, X, H, t1, T2, r, b, v) - S * Exp((b - r) * T2) * z6 + X * Exp(-r * T2) * z2
    ElseIf TypeFlag = "poB1" Then  '// put out type B1
        PartialTimeBarrier = PartialTimeBarrier("coB1", S, X, H, t1, T2, r, b, v) - S * Exp((b - r) * T2) * z8 + X * Exp(-r * T2) * z4 - S * Exp((b - r) * T2) * z7 + X * Exp(-r * T2) * z3
    ElseIf TypeFlag = "pdoB2" Then  '// put down-and-out type B2
        PartialTimeBarrier = PartialTimeBarrier("cdoB2", S, X, H, t1, T2, r, b, v) - S * Exp((b - r) * T2) * z7 + X * Exp(-r * T2) * z3
    ElseIf TypeFlag = "puoB2" Then  '// put up-and-out type B2
        PartialTimeBarrier = PartialTimeBarrier("cuoB2", S, X, H, t1, T2, r, b, v) - S * Exp((b - r) * T2) * z8 + X * Exp(-r * T2) * z4
    End If
End Function


'// Two asset barrier options
Public Function TwoAssetBarrier(TypeFlag As String, S1 As Double, S2 As Double, X As Double, H As Double, _
                T As Double, r As Double, b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double
    
    Dim d1 As Double, d2 As Double, d3 As Double, d4 As Double
    Dim e1 As Double, e2 As Double, e3 As Double, e4 As Double
    Dim mu1 As Double, mu2 As Double
    Dim eta As Integer    'Binary variable: 1 for call options and -1 for put options
    Dim phi As Integer    'Binary variable: 1 for up options and -1 for down options
    Dim KnockOutValue As Double
    
    mu1 = b1 - v1 ^ 2 / 2
    mu2 = b2 - v2 ^ 2 / 2
    
    d1 = (Log(S1 / X) + (mu1 + v1 ^ 2 / 2) * T) / (v1 * Sqr(T))
    d2 = d1 - v1 * Sqr(T)
    d3 = d1 + 2 * rho * Log(H / S2) / (v2 * Sqr(T))
    d4 = d2 + 2 * rho * Log(H / S2) / (v2 * Sqr(T))
    e1 = (Log(H / S2) - (mu2 + rho * v1 * v2) * T) / (v2 * Sqr(T))
    e2 = e1 + rho * v1 * Sqr(T)
    e3 = e1 - 2 * Log(H / S2) / (v2 * Sqr(T))
    e4 = e2 - 2 * Log(H / S2) / (v2 * Sqr(T))
   
    If TypeFlag = "cuo" Or TypeFlag = "cui" Then
        eta = 1: phi = 1
    ElseIf TypeFlag = "cdo" Or TypeFlag = "cdi" Then
        eta = 1: phi = -1
    ElseIf TypeFlag = "puo" Or TypeFlag = "pui" Then
        eta = -1: phi = 1
    ElseIf TypeFlag = "pdo" Or TypeFlag = "pdi" Then
        eta = -1: phi = -1
    End If
    KnockOutValue = eta * S1 * Exp((b1 - r) * T) * (CBND(eta * d1, phi * e1, -eta * phi * rho) _
    - Exp(2 * (mu2 + rho * v1 * v2) * Log(H / S2) / v2 ^ 2) * CBND(eta * d3, phi * e3, -eta * phi * rho)) - eta * Exp(-r * T) * X * (CBND(eta * d2, phi * e2, -eta * phi * rho) _
    - Exp(2 * mu2 * Log(H / S2) / v2 ^ 2) * CBND(eta * d4, phi * e4, -eta * phi * rho))
    If TypeFlag = "cuo" Or TypeFlag = "cdo" Or TypeFlag = "puo" Or TypeFlag = "pdo" Then
        TwoAssetBarrier = KnockOutValue
    ElseIf TypeFlag = "cui" Or TypeFlag = "cdi" Then
        TwoAssetBarrier = GBlackScholes("c", S1, X, T, r, b1, v1) - KnockOutValue
    ElseIf TypeFlag = "pui" Or TypeFlag = "pdi" Then
        TwoAssetBarrier = GBlackScholes("p", S1, X, T, r, b1, v1) - KnockOutValue
    End If
    
End Function


'// Partial-time two asset barrier options
Public Function PartialTimeTwoAssetBarrier(TypeFlag As String, S1 As Double, S2 As Double, X As Double, H As Double, t1 As Double, T2 As Double, _
                r As Double, b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double

    Dim d1 As Double, d2 As Double
    Dim d3 As Double, d4 As Double
    Dim e1 As Double, e2 As Double
    Dim e3 As Double, e4 As Double
    Dim mu1 As Double, mu2 As Double
    Dim OutBarrierValue As Double

    Dim eta As Integer
    Dim phi As Integer

    If TypeFlag = "cdo" Or TypeFlag = "pdo" Or TypeFlag = "cdi" Or TypeFlag = "pdi" Then
        phi = -1
    Else
        phi = 1
    End If
    
    If TypeFlag = "cdo" Or TypeFlag = "cuo" Or TypeFlag = "cdi" Or TypeFlag = "cui" Then
        eta = 1
    Else
        eta = -1
    End If
    mu1 = b1 - v1 ^ 2 / 2
    mu2 = b2 - v2 ^ 2 / 2
    d1 = (Log(S1 / X) + (mu1 + v1 ^ 2) * T2) / (v1 * Sqr(T2))
    d2 = d1 - v1 * Sqr(T2)
    d3 = d1 + 2 * rho * Log(H / S2) / (v2 * Sqr(T2))
    d4 = d2 + 2 * rho * Log(H / S2) / (v2 * Sqr(T2))
    e1 = (Log(H / S2) - (mu2 + rho * v1 * v2) * t1) / (v2 * Sqr(t1))
    e2 = e1 + rho * v1 * Sqr(t1)
    e3 = e1 - 2 * Log(H / S2) / (v2 * Sqr(t1))
    e4 = e2 - 2 * Log(H / S2) / (v2 * Sqr(t1))

    OutBarrierValue = eta * S1 * Exp((b1 - r) * T2) * (CBND(eta * d1, phi * e1, -eta * phi * rho * Sqr(t1 / T2)) - Exp(2 * Log(H / S2) * (mu2 + rho * v1 * v2) / (v2 ^ 2)) _
    * CBND(eta * d3, phi * e3, -eta * phi * rho * Sqr(t1 / T2))) _
    - eta * Exp(-r * T2) * X * (CBND(eta * d2, phi * e2, -eta * phi * rho * Sqr(t1 / T2)) - Exp(2 * Log(H / S2) * mu2 / (v2 ^ 2)) _
    * CBND(eta * d4, phi * e4, -eta * phi * rho * Sqr(t1 / T2)))

    If TypeFlag = "cdo" Or TypeFlag = "cuo" Or TypeFlag = "pdo" Or TypeFlag = "puo" Then
        PartialTimeTwoAssetBarrier = OutBarrierValue
     ElseIf TypeFlag = "cui" Or TypeFlag = "cdi" Then
        PartialTimeTwoAssetBarrier = GBlackScholes("c", S1, X, T2, r, b1, v1) - OutBarrierValue
    ElseIf TypeFlag = "pui" Or TypeFlag = "pdi" Then
        PartialTimeTwoAssetBarrier = GBlackScholes("p", S1, X, T2, r, b1, v1) - OutBarrierValue
    End If
End Function



'// Look-barrier options
Public Function LookBarrier(TypeFlag As String, S As Double, X As Double, H As Double, t1 As Double, T2 As Double, r As Double, b As Double, v As Double) As Double

    Dim hh As Double
    Dim K As Double, mu1 As Double, mu2 As Double
    Dim rho As Double, eta As Double, m As Double
    Dim g1 As Double, g2 As Double
    Dim OutValue As Double, part1 As Double, part2 As Double, part3 As Double, part4 As Double
    
    hh = Log(H / S)
    K = Log(X / S)
    mu1 = b - v ^ 2 / 2
    mu2 = b + v ^ 2 / 2
    rho = Sqr(t1 / T2)
    
    If TypeFlag = "cuo" Or TypeFlag = "cui" Then
        eta = 1
        m = Min(hh, K)
    ElseIf TypeFlag = "pdo" Or TypeFlag = "pdi" Then
        eta = -1
        m = Max(hh, K)
    End If
    
    g1 = (CND(eta * (hh - mu2 * t1) / (v * Sqr(t1))) - Exp(2 * mu2 * hh / v ^ 2) * CND(eta * (-hh - mu2 * t1) / (v * Sqr(t1)))) _
        - (CND(eta * (m - mu2 * t1) / (v * Sqr(t1))) - Exp(2 * mu2 * hh / v ^ 2) * CND(eta * (m - 2 * hh - mu2 * t1) / (v * Sqr(t1))))
    g2 = (CND(eta * (hh - mu1 * t1) / (v * Sqr(t1))) - Exp(2 * mu1 * hh / v ^ 2) * CND(eta * (-hh - mu1 * t1) / (v * Sqr(t1)))) _
        - (CND(eta * (m - mu1 * t1) / (v * Sqr(t1))) - Exp(2 * mu1 * hh / v ^ 2) * CND(eta * (m - 2 * hh - mu1 * t1) / (v * Sqr(t1))))

    part1 = S * Exp((b - r) * T2) * (1 + v ^ 2 / (2 * b)) * (CBND(eta * (m - mu2 * t1) / (v * Sqr(t1)), eta * (-K + mu2 * T2) / (v * Sqr(T2)), -rho) - Exp(2 * mu2 * hh / v ^ 2) _
        * CBND(eta * (m - 2 * hh - mu2 * t1) / (v * Sqr(t1)), eta * (2 * hh - K + mu2 * T2) / (v * Sqr(T2)), -rho))
    part2 = -Exp(-r * T2) * X * (CBND(eta * (m - mu1 * t1) / (v * Sqr(t1)), eta * (-K + mu1 * T2) / (v * Sqr(T2)), -rho) _
        - Exp(2 * mu1 * hh / v ^ 2) * CBND(eta * (m - 2 * hh - mu1 * t1) / (v * Sqr(t1)), eta * (2 * hh - K + mu1 * T2) / (v * Sqr(T2)), -rho))
    part3 = -Exp(-r * T2) * v ^ 2 / (2 * b) * (S * (S / X) ^ (-2 * b / v ^ 2) * CBND(eta * (m + mu1 * t1) / (v * Sqr(t1)), eta * (-K - mu1 * T2) / (v * Sqr(T2)), -rho) _
        - H * (H / X) ^ (-2 * b / v ^ 2) * CBND(eta * (m - 2 * hh + mu1 * t1) / (v * Sqr(t1)), eta * (2 * hh - K - mu1 * T2) / (v * Sqr(T2)), -rho))
    part4 = S * Exp((b - r) * T2) * ((1 + v ^ 2 / (2 * b)) * CND(eta * mu2 * (T2 - t1) / (v * Sqr(T2 - t1))) + Exp(-b * (T2 - t1)) * (1 - v ^ 2 / (2 * b)) _
        * CND(eta * (-mu1 * (T2 - t1)) / (v * Sqr(T2 - t1)))) * g1 - Exp(-r * T2) * X * g2
    OutValue = eta * (part1 + part2 + part3 + part4)

    If TypeFlag = "cuo" Or TypeFlag = "pdo" Then
        LookBarrier = OutValue
    ElseIf TypeFlag = "cui" Then
        LookBarrier = PartialFixedLB("c", S, X, t1, T2, r, b, v) - OutValue
    ElseIf TypeFlag = "pdi" Then
        LookBarrier = PartialFixedLB("p", S, X, t1, T2, r, b, v) - OutValue
    End If
End Function


'// Discrete barrier monitoring adjustment
Public Function DiscreteAdjustedBarrier(S As Double, H As Double, v As Double, dt As Double) As Double

    If H > S Then
        DiscreteAdjustedBarrier = H * Exp(0.5826 * v * Sqr(dt))
    ElseIf H < S Then
        DiscreteAdjustedBarrier = H * Exp(-0.5826 * v * Sqr(dt))
    End If
End Function


'// Soft barrier options
Public Function SoftBarrier(TypeFlag As String, S As Double, X As Double, _
                            L As Double, U As Double, T As Double, r As Double, b As Double, v As Double) As Double

    Dim mu As Double
    Dim d1 As Double, d2 As Double
    Dim d3 As Double, d4 As Double
    Dim e1 As Double, e2 As Double
    Dim e3 As Double, e4 As Double
    Dim lambda1 As Double, lambda2 As Double
    Dim Value As Double, eta As Integer
    
    If TypeFlag = "cdi" Or TypeFlag = "cdo" Then
        eta = 1
    Else
        eta = -1
    End If
    
    mu = (b + v ^ 2 / 2) / v ^ 2
    lambda1 = Exp(-1 / 2 * v ^ 2 * T * (mu + 0.5) * (mu - 0.5))
    lambda2 = Exp(-1 / 2 * v ^ 2 * T * (mu - 0.5) * (mu - 1.5))
    d1 = Log(U ^ 2 / (S * X)) / (v * Sqr(T)) + mu * v * Sqr(T)
    d2 = d1 - (mu + 0.5) * v * Sqr(T)
    d3 = Log(U ^ 2 / (S * X)) / (v * Sqr(T)) + (mu - 1) * v * Sqr(T)
    d4 = d3 - (mu - 0.5) * v * Sqr(T)
    e1 = Log(L ^ 2 / (S * X)) / (v * Sqr(T)) + mu * v * Sqr(T)
    e2 = e1 - (mu + 0.5) * v * Sqr(T)
    e3 = Log(L ^ 2 / (S * X)) / (v * Sqr(T)) + (mu - 1) * v * Sqr(T)
    e4 = e3 - (mu - 0.5) * v * Sqr(T)
    
    Value = eta * 1 / (U - L) * (S * Exp((b - r) * T) * S ^ (-2 * mu) _
    * (S * X) ^ (mu + 0.5) / (2 * (mu + 0.5)) _
    * ((U ^ 2 / (S * X)) ^ (mu + 0.5) * CND(eta * d1) - lambda1 * CND(eta * d2) _
    - (L ^ 2 / (S * X)) ^ (mu + 0.5) * CND(eta * e1) + lambda1 * CND(eta * e2)) _
    - X * Exp(-r * T) * S ^ (-2 * (mu - 1)) _
    * (S * X) ^ (mu - 0.5) / (2 * (mu - 0.5)) _
    * ((U ^ 2 / (S * X)) ^ (mu - 0.5) * CND(eta * d3) - lambda2 * CND(eta * d4) _
    - (L ^ 2 / (S * X)) ^ (mu - 0.5) * CND(eta * e3) + lambda2 * CND(eta * e4)))
    
    If TypeFlag = "cdi" Or TypeFlag = "pui" Then
        SoftBarrier = Value
    ElseIf TypeFlag = "cdo" Then
        SoftBarrier = GBlackScholes("c", S, X, T, r, b, v) - Value
    ElseIf TypeFlag = "puo" Then
        SoftBarrier = GBlackScholes("p", S, X, T, r, b, v) - Value
    End If
End Function


'// Gap options
Public Function GapOption(CallPutFlag As String, S As Double, X1 As Double, X2 As Double, _
                T As Double, r As Double, b As Double, v As Double) As Double

    Dim d1 As Double, d2 As Double

    d1 = (Log(S / X1) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
    
    If CallPutFlag = "c" Then
        GapOption = S * Exp((b - r) * T) * CND(d1) - X2 * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        GapOption = X2 * Exp(-r * T) * CND(-d2) - S * Exp((b - r) * T) * CND(-d1)
    End If
End Function


'// Cash-or-nothing options
Public Function CashOrNothing(CallPutFlag As String, S As Double, X As Double, K As Double, T As Double, _
                r As Double, b As Double, v As Double) As Double

    Dim d As Double

    d = (Log(S / X) + (b - v ^ 2 / 2) * T) / (v * Sqr(T))

    If CallPutFlag = "c" Then
        CashOrNothing = K * Exp(-r * T) * CND(d)
    ElseIf CallPutFlag = "p" Then
        CashOrNothing = K * Exp(-r * T) * CND(-d)
    End If
End Function


'// Two asset cash-or-nothing options
Public Function TwoAssetCashOrNothing(TypeFlag As Integer, S1 As Double, S2 As Double, X1 As Double, X2 As Double, K As Double, T As Double, r As Double, _
                b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double) As Double
    
    Dim d1 As Double, d2 As Double
                                   
    d1 = (Log(S1 / X1) + (b1 - v1 ^ 2 / 2) * T) / (v1 * Sqr(T))
    d2 = (Log(S2 / X2) + (b2 - v2 ^ 2 / 2) * T) / (v2 * Sqr(T))
                                
    If TypeFlag = 1 Then
        TwoAssetCashOrNothing = K * Exp(-r * T) * CBND(d1, d2, rho)
    ElseIf TypeFlag = 2 Then
        TwoAssetCashOrNothing = K * Exp(-r * T) * CBND(-d1, -d2, rho)
    ElseIf TypeFlag = 3 Then
        TwoAssetCashOrNothing = K * Exp(-r * T) * CBND(d1, -d2, -rho)
    ElseIf TypeFlag = 4 Then
        TwoAssetCashOrNothing = K * Exp(-r * T) * CBND(-d1, d2, -rho)
    End If
End Function


'// Asset-or-nothing options
Public Function AssetOrNothing(CallPutFlag As String, S As Double, X As Double, T As Double, r As Double, _
                b As Double, v As Double) As Double

    Dim d As Double
    
    d = (Log(S / X) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    
    If CallPutFlag = "c" Then
        AssetOrNothing = S * Exp((b - r) * T) * CND(d)
    ElseIf CallPutFlag = "p" Then
        AssetOrNothing = S * Exp((b - r) * T) * CND(-d)
    End If
End Function


'// Supershare options
Public Function SuperShare(S As Double, XL As Double, XH As Double, T As Double, _
                r As Double, b As Double, v As Double) As Double
 
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(S / XL) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = (Log(S / XH) + (b + v ^ 2 / 2) * T) / (v * Sqr(T))

    SuperShare = S * Exp((b - r) * T) / XL * (CND(d1) - CND(d2))
End Function


'// Binary barrier options
Public Function BinaryBarrier(TypeFlag As Integer, S As Double, X As Double, H As Double, K As Double, _
                T As Double, r As Double, b As Double, v As Double, eta As Integer, phi As Integer) As Double

    '// TypeFlag:  Value 1 to 28 dependent on binary option type,
    '//            look in the book for spesifications.
    
    Dim X1 As Double, X2 As Double
    Dim y1 As Double, y2 As Double
    Dim Z As Double, mu As Double, lambda As Double
    Dim a1 As Double, a2 As Double, a3 As Double, a4 As Double, a5 As Double
    Dim b1 As Double, b2 As Double, b3 As Double, b4 As Double

    mu = (b - v ^ 2 / 2) / v ^ 2
    lambda = Sqr(mu ^ 2 + 2 * r / v ^ 2)
    X1 = Log(S / X) / (v * Sqr(T)) + (mu + 1) * v * Sqr(T)
    X2 = Log(S / H) / (v * Sqr(T)) + (mu + 1) * v * Sqr(T)
    y1 = Log(H ^ 2 / (S * X)) / (v * Sqr(T)) + (mu + 1) * v * Sqr(T)
    y2 = Log(H / S) / (v * Sqr(T)) + (mu + 1) * v * Sqr(T)
    Z = Log(H / S) / (v * Sqr(T)) + lambda * v * Sqr(T)
    
    a1 = S * Exp((b - r) * T) * CND(phi * X1)
    b1 = K * Exp(-r * T) * CND(phi * X1 - phi * v * Sqr(T))
    a2 = S * Exp((b - r) * T) * CND(phi * X2)
    b2 = K * Exp(-r * T) * CND(phi * X2 - phi * v * Sqr(T))
    a3 = S * Exp((b - r) * T) * (H / S) ^ (2 * (mu + 1)) * CND(eta * y1)
    b3 = K * Exp(-r * T) * (H / S) ^ (2 * mu) * CND(eta * y1 - eta * v * Sqr(T))
    a4 = S * Exp((b - r) * T) * (H / S) ^ (2 * (mu + 1)) * CND(eta * y2)
    b4 = K * Exp(-r * T) * (H / S) ^ (2 * mu) * CND(eta * y2 - eta * v * Sqr(T))
    a5 = K * ((H / S) ^ (mu + lambda) * CND(eta * Z) + (H / S) ^ (mu - lambda) * CND(eta * Z - 2 * eta * lambda * v * Sqr(T)))
    
    If X > H Then
        Select Case TypeFlag
            Case Is < 5
                BinaryBarrier = a5
            Case Is < 7
                BinaryBarrier = b2 + b4
            Case Is < 9
                BinaryBarrier = a2 + a4
            Case Is < 11
                BinaryBarrier = b2 - b4
            Case Is < 13
                BinaryBarrier = a2 - a4
            Case Is = 13
                BinaryBarrier = b3
            Case Is = 14
                BinaryBarrier = b3
            Case Is = 15
                BinaryBarrier = a3
            Case Is = 16
                BinaryBarrier = a1
            Case Is = 17
                BinaryBarrier = b2 - b3 + b4
            Case Is = 18
                BinaryBarrier = b1 - b2 + b4
            Case Is = 19
                BinaryBarrier = a2 - a3 + a4
            Case Is = 20
                BinaryBarrier = a1 - a2 + a3
            Case Is = 21
                BinaryBarrier = b1 - b3
            Case Is = 22
                BinaryBarrier = 0
            Case Is = 23
                BinaryBarrier = a1 - a3
            Case Is = 24
               BinaryBarrier = 0
            Case Is = 25
                BinaryBarrier = b1 - b2 + b3 - b4
            Case Is = 26
                BinaryBarrier = b2 - b4
            Case Is = 27
                BinaryBarrier = a1 - a2 + a3 - a4
            Case Is = 28
                BinaryBarrier = a2 - a4
        End Select
    ElseIf X < H Then
        Select Case TypeFlag
            Case Is < 5
                BinaryBarrier = a5
            Case Is < 7
                BinaryBarrier = b2 + b4
            Case Is < 9
                BinaryBarrier = a2 + a4
            Case Is < 11
                BinaryBarrier = b2 - b4
            Case Is < 13
                BinaryBarrier = a2 - a4
            Case Is = 13
                BinaryBarrier = b1 - b2 + b4
            Case Is = 14
                BinaryBarrier = b2 - b3 + b4
            Case Is = 15
                BinaryBarrier = a1 - a2 + a4
            Case Is = 16
                BinaryBarrier = a2 - a3 + a4
            Case Is = 17
                BinaryBarrier = b1
            Case Is = 18
                BinaryBarrier = b3
            Case Is = 19
                BinaryBarrier = a1
            Case Is = 20
                BinaryBarrier = a3
            Case Is = 21
                BinaryBarrier = b2 - b4
            Case Is = 22
                BinaryBarrier = b1 - b2 + b3 - b4
            Case Is = 23
                BinaryBarrier = a2 - a4
            Case Is = 24
                BinaryBarrier = a1 - a2 + a3 - a4
            Case Is = 25
                BinaryBarrier = 0
            Case Is = 26
                BinaryBarrier = b1 - b3
            Case Is = 27
                BinaryBarrier = 0
            Case Is = 28
                BinaryBarrier = a1 - a3
        End Select
    End If
End Function


'// Geometric average rate option
Public Function GeometricAverageRateOption(CallPutFlag As String, S As Double, SA As Double, X As Double, _
                T As Double, T2 As Double, r As Double, b As Double, v As Double) As Double


    Dim t1 As Double 'Observed or realized time period
    Dim bA As Double, vA As Double
    
    bA = 1 / 2 * (b - v ^ 2 / 6)
    vA = v / Sqr(3)

    t1 = T - T2
    
    If t1 > 0 Then
        X = (t1 + T2) / T2 * X - t1 / T2 * SA
        GeometricAverageRateOption = GBlackScholes(CallPutFlag, S, X, T2, r, bA, vA) * T2 / (t1 + T2)
    ElseIf t1 = 0 Then
        GeometricAverageRateOption = GBlackScholes(CallPutFlag, S, X, T, r, bA, vA)
    End If

End Function


'// Arithmetic average rate option
Public Function TurnbullWakemanAsian(CallPutFlag As String, S As Double, SA As Double, X As Double, _
            T As Double, T2 As Double, tau As Double, r As Double, b As Double, v As Double) As Double

    Dim m1 As Double, m2 As Double, t1 As Double
    Dim bA As Double, vA As Double
    
    m1 = (Exp(b * T) - Exp(b * tau)) / (b * (T - tau))
    m2 = 2 * Exp((2 * b + v ^ 2) * T) / ((b + v ^ 2) * (2 * b + v ^ 2) * (T - tau) ^ 2) _
    + 2 * Exp((2 * b + v ^ 2) * tau) / (b * (T - tau) ^ 2) * (1 / (2 * b + v ^ 2) - Exp(b * (T - tau)) / (b + v ^ 2))
    
    bA = Log(m1) / T
    vA = Sqr(Log(m2) / T - 2 * bA)
    t1 = T - T2
    
    If t1 > 0 Then
        X = T / T2 * X - t1 / T2 * SA
        TurnbullWakemanAsian = GBlackScholes(CallPutFlag, S, X, T2, r, bA, vA) * T2 / T
    Else
        TurnbullWakemanAsian = GBlackScholes(CallPutFlag, S, X, T2, r, bA, vA)
    End If
End Function


'// Arithmetic average rate option
Public Function LevyAsian(CallPutFlag As String, S As Double, SA As Double, X As Double, _
                T As Double, T2 As Double, r As Double, b As Double, v As Double) As Double

    Dim SE As Double
    Dim m As Double, d As Double
    Dim Sv As Double, XStar As Double
    Dim d1 As Double, d2 As Double

    SE = S / (T * b) * (Exp((b - r) * T2) - Exp(-r * T2))
    m = 2 * S ^ 2 / (b + v ^ 2) * ((Exp((2 * b + v ^ 2) * T2) - 1) / (2 * b + v ^ 2) - (Exp(b * T2) - 1) / b)
    d = m / (T ^ 2)
    Sv = Log(d) - 2 * (r * T2 + Log(SE))
    XStar = X - (T - T2) / T * SA
    d1 = 1 / Sqr(Sv) * (Log(d) / 2 - Log(XStar))
    d2 = d1 - Sqr(Sv)
    
    If CallPutFlag = "c" Then
        LevyAsian = SE * CND(d1) - XStar * Exp(-r * T2) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        LevyAsian = (SE * CND(d1) - XStar * Exp(-r * T2) * CND(d2)) - SE + XStar * Exp(-r * T2)
    End If
End Function


'// Foreign equity option struck in domestic currency
Public Function ForEquOptInDomCur(CallPutFlag As String, E As Double, S As Double, X As Double, T As Double, _
    r As Double, q As Double, vS As Double, vE As Double, rho As Double) As Double

    Dim v As Double, d1 As Double, d2 As Double

    v = Sqr(vE ^ 2 + vS ^ 2 + 2 * rho * vE * vS)
    d1 = (Log(E * S / X) + (r - q + v ^ 2 / 2) * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)
   
    If CallPutFlag = "c" Then
        ForEquOptInDomCur = E * S * Exp(-q * T) * CND(d1) - X * Exp(-r * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        ForEquOptInDomCur = X * Exp(-r * T) * CND(-d2) - E * S * Exp(-q * T) * CND(-d1)
    End If
End Function


'// Fixed exchange rate foreign equity options-- Quantos
Public Function Quanto(CallPutFlag As String, Ep As Double, S As Double, X As Double, T As Double, r As Double, _
                rf As Double, q As Double, vS As Double, vE As Double, rho As Double) As Double
    
    Dim d1 As Double, d2 As Double

    d1 = (Log(S / X) + (rf - q - rho * vS * vE + vS ^ 2 / 2) * T) / (vS * Sqr(T))
    d2 = d1 - vS * Sqr(T)
   
    If CallPutFlag = "c" Then
        Quanto = Ep * (S * Exp((rf - r - q - rho * vS * vE) * T) * CND(d1) - X * Exp(-r * T) * CND(d2))
    ElseIf CallPutFlag = "p" Then
        Quanto = Ep * (X * Exp(-r * T) * CND(-d2) - S * Exp((rf - r - q - rho * vS * vE) * T) * CND(-d1))
    End If
End Function


'// Equity linked foreign exchange option
Public Function EquityLinkedFXO(CallPutFlag As String, E As Double, S As Double, X As Double, T As Double, r As Double, _
                rf As Double, q As Double, vS As Double, vE As Double, rho As Double) As Double

    Dim d1 As Double, d2 As Double
    
    d1 = (Log(E / X) + (r - rf + rho * vS * vE + vE ^ 2 / 2) * T) / (vE * Sqr(T))
    d2 = d1 - vE * Sqr(T)
     
    If CallPutFlag = "c" Then
        EquityLinkedFXO = E * S * Exp(-q * T) * CND(d1) - X * S * Exp((rf - r - q - rho * vS * vE) * T) * CND(d2)
    ElseIf CallPutFlag = "p" Then
        EquityLinkedFXO = X * S * Exp((rf - r - q - rho * vS * vE) * T) * CND(-d2) - E * S * Exp(-q * T) * CND(-d1)
    End If
End Function


'// Takeover foreign exchange options
Public Function TakeoverFXoption(v As Double, b As Double, E As Double, X As Double, T As Double, r As Double, rf As Double, vV As Double, vE As Double, rho As Double) As Double
    
    Dim a1 As Double, a2 As Double
    
    a1 = (Log(v / b) + (rf - rho * vE * vV - vV ^ 2 / 2) * T) / (vV * Sqr(T))
    a2 = (Log(E / X) + (r - rf - vE ^ 2 / 2) * T) / (vE * Sqr(T))
    
    TakeoverFXoption = b * (E * Exp(-rf * T) * CBND(a2 + vE * Sqr(T), -a1 - rho * vE * Sqr(T), -rho) _
    - X * Exp(-r * T) * CBND(-a1, a2, -rho))

End Function


'//  Black-76 European swaption
Public Function Swaption(CallPutFlag As String, t1 As Double, m As Double, F As Double, X As Double, T As Double, _
                r As Double, v As Double) As Double
 
    Dim d1 As Double, d2 As Double
    
    d1 = (Log(F / X) + v ^ 2 / 2 * T) / (v * Sqr(T))
    d2 = d1 - v * Sqr(T)

    If CallPutFlag = "c" Then 'Payer swaption
        Swaption = ((1 - 1 / (1 + F / m) ^ (t1 * m)) / F) * Exp(-r * T) * (F * CND(d1) - X * CND(d2))
    ElseIf CallPutFlag = "p" Then  'Receiver swaption
        Swaption = ((1 - 1 / (1 + F / m) ^ (t1 * m)) / F) * Exp(-r * T) * (X * CND(-d2) - F * CND(-d1))
    End If

End Function


'// Vasicek: options on zero coupon bonds
Function VasicekBondOption(CallPutFlag As String, F As Double, X As Double, tau As Double, T As Double, _
        r As Double, theta As Double, kappa As Double, v As Double) As Double

  
    Dim PtT As Double, Pt_tau As Double
    Dim H As Double, vp As Double

    X = X / F
    PtT = VasicekBondPrice(0, T, r, theta, kappa, v)
    Pt_tau = VasicekBondPrice(0, tau, r, theta, kappa, v)
    vp = Sqr(v ^ 2 * (1 - Exp(-2 * kappa * T)) / (2 * kappa)) * (1 - Exp(-kappa * (tau - T))) / kappa
   
    H = 1 / vp * Log(Pt_tau / (PtT * X)) + vp / 2
    
    If CallPutFlag = "c" Then
        VasicekBondOption = F * (Pt_tau * CND(H) - X * PtT * CND(H - vp))
    Else
        VasicekBondOption = F * (X * PtT * CND(-H + vp) - Pt_tau * CND(-H))
    End If
End Function


'// Vasicek: value zero coupon bond
Public Function VasicekBondPrice(t1 As Double, T As Double, r As Double, theta As Double, kappa As Double, v As Double) As Double
    Dim BtT As Double, AtT As Double, PtT As Double

    BtT = (1 - Exp(-kappa * (T - t1))) / kappa
    AtT = Exp((BtT - T + t1) * (kappa ^ 2 * theta - v ^ 2 / 2) / kappa ^ 2 - v ^ 2 * BtT ^ 2 / (4 * kappa))
    PtT = AtT * Exp(-BtT * r)
    VasicekBondPrice = PtT
End Function




Option Explicit     'Requires that all variables to be declared explicitly.
Option Compare Text 'Uppercase letters to be equivalent to lowercase letters.

Option Base 0       'The "Option Base" statement allows to specify 0 or 1 as the
                    'default first index of arrays.


                    
Public Function Max(X, y)
            Max = Application.Max(X, y)
End Function

Public Function Min(X, y)
            Min = Application.Min(X, y)
End Function

'// Cox-Ross-Rubinstein binomial tree
Public Function CRRBinomial(AmeEurFlag As String, CallPutFlag As String, S As Double, X As Double, T As Double, _
                r As Double, b As Double, v As Double, n As Integer) As Double
                

    Dim OptionValue() As Double
    Dim u As Double, d As Double, p As Double
    Dim dt As Double, Df As Double
    Dim i As Integer, j As Integer, z As Integer
    
    ReDim OptionValue(n + 1)
    
    If CallPutFlag = "c" Then
        z = 1
        ElseIf CallPutFlag = "p" Then
        z = -1
    End If
    
    dt = T / n
    u = Exp(v * Sqr(dt))
    d = 1 / u
    p = (Exp(b * dt) - d) / (u - d)
    Df = Exp(-r * dt)
    
    For i = 0 To n
         OptionValue(i) = Max(0, z * (S * u ^ i * d ^ (n - i) - X))
    Next
    
    For j = n - 1 To 0 Step -1:
        For i = 0 To j
            If AmeEurFlag = "e" Then
                OptionValue(i) = (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df
            ElseIf AmeEurFlag = "a" Then
                OptionValue(i) = Max((z * (S * u ^ i * d ^ (Abs(i - j)) - X)), _
                (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df)
            End If
        Next
    Next
    CRRBinomial = OptionValue(0)
End Function


'// European and American barrier options in binomial trees
Public Function BarrierBinomial(AmeEurFlag As String, TypeFlag As String, S As Double, X As Double, H As Double, T As Double, _
                r As Double, b As Double, v As Double, n As Integer) As Double

    Dim AssetPrice As Double
    Dim OptionValue() As Double
    Dim dt As Double, Df As Double
    Dim u As Double, d As Double, p As Double
    Dim i As Integer, j As Integer
    Dim phi As Integer, z As Integer
    
     For i = 1 To 100
        If n < Int((i ^ 2 * v ^ 2 * T) / (Log(S / H)) ^ 2) Then
            n = Int((i ^ 2 * v ^ 2 * T) / (Log(S / H)) ^ 2)
            Exit For
        End If
    Next
    
    ReDim OptionValue(n + 1)
    
    If TypeFlag = "cuo" Or TypeFlag = "cdo" Then
        z = 1
        ElseIf TypeFlag = "puo" Or TypeFlag = "pdo" Then
        z = -1
    End If
    If TypeFlag = "cuo" Or TypeFlag = "puo" Then
        phi = 1
        ElseIf TypeFlag = "cdo" Or TypeFlag = "pdo" Then
        phi = -1
    End If
    
    dt = T / n
    u = Exp(v * Sqr(dt))
    d = 1 / u
    p = (Exp(b * dt) - d) / (u - d)
    Df = Exp(-r * dt)
    
    For i = 0 To n
         OptionValue(i) = Max(0, z * (S * u ^ i * d ^ (n - i) - X))
    Next
    
    For j = n - 1 To 0 Step -1:
        For i = 0 To j
            AssetPrice = S * u ^ i * d ^ Abs(i - j)
            
            If phi = 1 Then
                If AmeEurFlag = "e" And AssetPrice < H Then
                    OptionValue(i) = (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df
                ElseIf AmeEurFlag = "a" And AssetPrice < H Then
                    OptionValue(i) = Max((z * (AssetPrice - X)), _
                    (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df)
                ElseIf AssetPrice >= H Then
                    OptionValue(i) = 0
                End If
            Else
                If AmeEurFlag = "e" And AssetPrice > H Then
                    OptionValue(i) = (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df
                ElseIf AmeEurFlag = "a" And AssetPrice > H Then
                    OptionValue(i) = Max((z * (AssetPrice - X)), _
                    (p * OptionValue(i + 1) + (1 - p) * OptionValue(i)) * Df)
                ElseIf AssetPrice <= H Then
                    OptionValue(i) = 0
                End If
            End If
        Next
    Next
    BarrierBinomial = OptionValue(0)
End Function


'// Convertible bonds in binomial trees
Public Function ConvertibleBond(AmeEurFlag As String, S As Double, X As Double, T2 As Double, t1 As Double, _
                r As Double, k As Double, q As Double, v As Double, F As Double, Coupon As Double, n As Integer) As Double
                
    'T2:            Time to maturity
    't1:            Time to to first conversion date, only used when the option is American
    'k              Issuer's credit spread above risk-free return (similar Treasuries).
    'q              Stock dividend yield
    'v:             Stock volatility
    'n:             Number of time-steps in the tree
    'F              Bond face value
    
    Dim OptionValue() As Double           'Option Value at each node
    Dim ConversionProbability() As Double 'Conversion probability at each node
    Dim dt As Double
    Dim u As Double, d As Double, p As Double
    Dim ru As Double, rd As Double
    Dim CouponValue As Double
    Dim CouponTime As Double    'Time to last coupon payment from today
    Dim Time As Double
    Dim StockPrice As Double
    Dim CRatio As Double    'Conversion Ratio= FaceValue/Strike
    Dim i As Integer, j As Integer, CouponPayment As Integer
    
    CRatio = F / X
    
    ReDim OptionValue(n + 1)
    ReDim ConversionProbability(n + 1)
    
    dt = T2 / n
    u = Exp(v * Sqr(dt))
    d = 1 / u
    p = (Exp((r - q) * dt) - d) / (u - d)
    
    For i = 0 To n
        StockPrice = S * u ^ i * d ^ (n - i)
        OptionValue(i) = Max(F + Coupon, CRatio * StockPrice)
        If OptionValue(i) = F + Coupon Then
            ConversionProbability(i) = 0
        Else
            ConversionProbability(i) = 1
        End If
    Next
    
    CouponTime = T2
    For j = n - 1 To 0 Step -1:
     
     If CouponTime <> T2 - Int((n - j) * dt) Then
        CouponPayment = 1
        CouponTime = T2 - Int((n - j) * dt)
        Time = CouponTime - dt * j
    Else
        CouponPayment = 0
     End If

        For i = 0 To j
            ru = ConversionProbability(i + 1) * r + (1 - ConversionProbability(i + 1)) * (r + k)
            rd = ConversionProbability(i) * r + (1 - ConversionProbability(i)) * (r + k)
            StockPrice = S * u ^ i * d ^ Abs(i - j)
            If CouponPayment = 1 Then
                CouponValue = Coupon * (p * Exp(-ru * Time) + (1 - p) * Exp(-rd * Time))
            Else
                CouponValue = 0
            End If
            
            If AmeEurFlag = "e" Or AmeEurFlag = "a" And dt * j <= t1 Then
                OptionValue(i) = CouponValue + p * OptionValue(i + 1) * Exp(-ru * dt) + (1 - p) * OptionValue(i) * Exp(-rd * dt)
                ConversionProbability(i) = p * ConversionProbability(i + 1) + (1 - p) * ConversionProbability(i)
            ElseIf AmeEurFlag = "a" Then
                OptionValue(i) = Max(CRatio * (StockPrice - X) + F, _
                CouponValue + p * OptionValue(i + 1) * Exp(-ru * dt) + (1 - p) * OptionValue(i) * Exp(-rd * dt))
                If OptionValue(i) = CRatio * (StockPrice - X) + F Then
                    ConversionProbability(i) = 1
                Else
                    ConversionProbability(i) = p * ConversionProbability(i + 1) + (1 - p) * ConversionProbability(i)
                End If
            End If
        Next
    Next
    ConvertibleBond = OptionValue(0)
End Function



'// Trinomial tree
Public Function TrinomialTree(AmeEurFlag As String, CallPutFlag As String, S As Double, X As Double, T As Double, _
                r As Double, b As Double, v As Double, n As Integer) As Double
                

    Dim OptionValue() As Double
    Dim dt As Double, u As Double, d As Double
    Dim pu As Double, pd As Double, pm As Double
    Dim i As Integer, j As Integer, z As Integer
    Dim Df As Double
    
    ReDim OptionValue(n * 2 + 1)
    
    If CallPutFlag = "c" Then
        z = 1
        ElseIf CallPutFlag = "p" Then
        z = -1
    End If
    
    dt = T / n
    u = Exp(v * Sqr(2 * dt))
    d = Exp(-v * Sqr(2 * dt))
    pu = ((Exp(b * dt / 2) - Exp(-v * Sqr(dt / 2))) / (Exp(v * Sqr(dt / 2)) - Exp(-v * Sqr(dt / 2)))) ^ 2
   
    pd = ((Exp(v * Sqr(dt / 2)) - Exp(b * dt / 2)) / (Exp(v * Sqr(dt / 2)) - Exp(-v * Sqr(dt / 2)))) ^ 2
    pm = 1 - pu - pd
    Df = Exp(-r * dt)
    
    For i = 0 To (2 * n)
         OptionValue(i) = Max(0, z * (S * u ^ Max(i - n, 0) * d ^ Max(n * 2 - n - i, 0) - X))
    Next
    For j = n - 1 To 0 Step -1
        For i = 0 To (j * 2)
            If AmeEurFlag = "e" Then
                OptionValue(i) = (pu * OptionValue(i + 2) + pm * OptionValue(i + 1) + pd * OptionValue(i)) * Df
            ElseIf AmeEurFlag = "a" Then
                OptionValue(i) = Max((z * (S * u ^ Max(i - j, 0) * d ^ Max(j * 2 - j - i, 0) - X)), _
                (pu * OptionValue(i + 2) + pm * OptionValue(i + 1) + pd * OptionValue(i)) * Df)
            End If
        Next
    Next
    TrinomialTree = OptionValue(0)
End Function


'// Three dimensional binomial tree
Public Function ThreeDimensionalBinomial(TypeFlag As Integer, AmeEurFlag As String, CallPutFlag As String, _
                S1 As Double, S2 As Double, Q1 As Double, Q2 As Double, X1 As Double, X2 As Double, T As Double, r As Double, b1 As Double, _
                b2 As Double, v1 As Double, v2 As Double, rho As Double, n As Integer) As Double
    
    Dim OptionValue() As Double
    Dim dt As Double, u As Double, d As Double
    Dim mu1 As Double, mu2 As Double
    Dim Y1 As Double, y2
    Dim NodeValueS1 As Double, NodeValueS2 As Double
    Dim i As Integer, j As Integer, m As Integer
    
    ReDim OptionValue(n + 1, n + 1)

    dt = T / n
    mu1 = b1 - v1 ^ 2 / 2
    mu2 = b2 - v2 ^ 2 / 2
    u = Exp(mu1 * dt + v1 * Sqr(dt))
    d = Exp(mu1 * dt - v1 * Sqr(dt))

    For j = 0 To n
        Y1 = (2 * j - n) * Sqr(dt)
        NodeValueS1 = S1 * u ^ j * d ^ (n - j)
        For i = 0 To n
            y2 = rho * Y1 + Sqr(1 - rho ^ 2) * (2 * i - n) * Sqr(dt)
            NodeValueS2 = S2 * Exp(mu2 * n * dt) * Exp(v2 * (rho * Y1 + Sqr(1 - rho ^ 2) * (2 * i - n) * Sqr(dt)))
            OptionValue(j, i) = PayoffFunction(TypeFlag, CallPutFlag, NodeValueS1, NodeValueS2, Q1, Q2, X1, X2)
        Next
    Next
    For m = n - 1 To 0 Step -1
        For j = 0 To m
            Y1 = (2 * j - m) * Sqr(dt)
            NodeValueS1 = S1 * u ^ j * d ^ (m - j)
            For i = 0 To m
                y2 = rho * Y1 + Sqr(1 - rho ^ 2) * (2 * i - m) * Sqr(dt)
                NodeValueS2 = S2 * Exp(mu2 * m * dt) * Exp(v2 * y2)
                OptionValue(j, i) = 1 / 4 * (OptionValue(j, i) + OptionValue(j + 1, i) + OptionValue(j, i + 1) _
                + OptionValue(j + 1, i + 1)) * Exp(-r * dt)
                If AmeEurFlag = "a" Then
                    OptionValue(j, i) = Max(OptionValue(j, i), PayoffFunction(TypeFlag, CallPutFlag, NodeValueS1, NodeValueS2, Q1, Q2, X1, X2))
                End If
            Next
        Next
    Next
    ThreeDimensionalBinomial = OptionValue(0, 0)
End Function
'// Payoff function used in three dimensional binomial tree
Public Function PayoffFunction(TypeFlag As Integer, CallPutFlag As String, S1 As Double, S2 As Double, _
        Q1 As Double, Q2 As Double, X1 As Double, X2 As Double) As Double

    Dim z As Integer
    
    If CallPutFlag = "c" Then
        z = 1
        ElseIf CallPutFlag = "p" Then
        z = -1
    End If

    If TypeFlag = 1 Then '// Spread option
        PayoffFunction = Max(0, z * (Q1 * S1 - Q2 * S2) - z * X1)
    ElseIf TypeFlag = 2 Then  '// Option on the maximum of two assets
        PayoffFunction = Max(0, z * Max(Q1 * S1, Q2 * S2) - z * X1)
    ElseIf TypeFlag = 3 Then  '// Option on the minimum of two assets
        PayoffFunction = Max(0, z * Min(Q1 * S1, Q2 * S2) - z * X1)
    ElseIf TypeFlag = 4 Then  '// Dual strike option
        PayoffFunction = Application.Max(0, z * (Q1 * S1 - X1), z * (Q2 * S2 - X2))
    ElseIf TypeFlag = 5 Then  '// Reverse-dual strike option
        PayoffFunction = Application.Max(0, z * (Q1 * S1 - X1), z * (X2 - Q2 * S2))
    ElseIf TypeFlag = 6 Then  '// Portfolio option
        PayoffFunction = Max(0, z * (Q1 * S1 + Q2 * S2) - z * X1)
    ElseIf TypeFlag = 7 Then  '// Exchange option
        PayoffFunction = Max(0, Max(0, Q1 * S1 - Q2 * S2))
    End If
End Function


'// Implied trinomial tree
Public Function ImpliedTrinomialTree(ReturnFlag As String, STEPn As Integer, STATEi As Integer, S As Double, X As Double, T As Double, _
                r As Double, b As Double, v As Double, Skew As Double, nSteps As Integer)
                
    Dim ArrowDebreu() As Double
    Dim LocalVolatility() As Double
    Dim UpProbability() As Double
    Dim DownProbability() As Double
    Dim OptionValueNode() As Double
    Dim dt As Double, u As Double, d As Double
    Dim Df As Double, pi As Double, qi As Double
    Dim Si1 As Double, Si As Double, Si2 As Double
    Dim vi As Double, Fj As Double, Fi As Double, Fo As Double
    Dim Sum As Double, OptionValue As Double
    Dim i As Integer, j As Integer, n As Integer, z As Integer
   
    ReDim OptionValueNode(nSteps * 2) As Double
    ReDim ArrowDebreu(nSteps, nSteps * 2) As Double
    ReDim UpProbability(nSteps - 1, nSteps * 2 - 2) As Double
    ReDim DownProbability(nSteps - 1, nSteps * 2 - 2) As Double
    ReDim LocalVolatility(nSteps - 1, nSteps * 2 - 2) As Double
   
    dt = T / nSteps
    u = Exp(v * Sqr(2 * dt))
    d = 1 / u
    Df = Exp(-r * dt)
    ArrowDebreu(0, 0) = 1
    For n = 0 To nSteps - 1
        For i = 0 To n * 2
            Sum = 0
            Si1 = S * u ^ Max(i - n, 0) * d ^ Max(n * 2 - n - i, 0)
            Si = Si1 * d
            Si2 = Si1 * u
            Fi = Si1 * Exp(b * dt)
            vi = v + (S - Si1) * Skew
            If i < (n * 2) / 2 + 1 Then
                For j = 0 To i - 1
                    Fj = S * u ^ Max(j - n, 0) * d ^ Max(n * 2 - n - j, 0) * Exp(b * dt)
                    Sum = Sum + ArrowDebreu(n, j) * (Si1 - Fj)
                Next
                OptionValue = TrinomialTree("e", "p", S, Si1, (n + 1) * dt, r, b, vi, n + 1)
                qi = (Exp(r * dt) * OptionValue - Sum) / (ArrowDebreu(n, i) * (Si1 - Si))
                pi = (Fi + qi * (Si1 - Si) - Si1) / (Si2 - Si1)
            Else
                OptionValue = TrinomialTree("e", "c", S, Si1, (n + 1) * dt, r, b, vi, n + 1)
                Sum = 0
                For j = i + 1 To n * 2
                    Fj = S * u ^ Max(j - n, 0) * d ^ Max(n * 2 - n - j, 0) * Exp(b * dt)
                    Sum = Sum + ArrowDebreu(n, j) * (Fj - Si1)
                Next
                pi = (Exp(r * dt) * OptionValue - Sum) / (ArrowDebreu(n, i) * (Si2 - Si1))
                qi = (Fi - pi * (Si2 - Si1) - Si1) / (Si - Si1)
            End If
            
            '// Replacing negative probabilities
            If pi < 0 Or pi > 1 Or qi < 0 Or qi > 1 Then
                If Fi > Si1 And Fi < Si2 Then
                    pi = 1 / 2 * ((Fi - Si1) / (Si2 - Si1) + (Fi - Si) / (Si2 - Si))
                    qi = 1 / 2 * ((Si2 - Fi) / (Si2 - Si))
                ElseIf Fi > Si And Fi < Si1 Then
                    pi = 1 / 2 * ((Fi - Si) / (Si2 - Si))
                    qi = 1 / 2 * ((Si2 - Fi) / (Si2 - Si) + (Si1 - Fi) / (Si1 - Si))
                End If
            End If
                DownProbability(n, i) = qi
                UpProbability(n, i) = pi
            '// Calculating local volatilities
                Fo = (pi * Si2 + qi * Si + (1 - pi - qi) * Si1) * Exp(b * dt)
                LocalVolatility(n, i) = Sqr((pi * (Si2 - Fo) ^ 2 + (1 - pi - qi) * _
                    (Si1 - Fo) ^ 2 + qi * (Si - Fo) ^ 2) / (Fo ^ 2 * dt))
                    
            '// Calculating Arrow-Debreu prices
            If n = 0 Then
                ArrowDebreu(n + 1, i) = qi * ArrowDebreu(n, i) * Df
                ArrowDebreu(n + 1, i + 1) = (1 - pi - qi) * ArrowDebreu(n, i) * Df
                ArrowDebreu(n + 1, i + 2) = pi * ArrowDebreu(n, i) * Df
            ElseIf n > 0 And i = 0 Then
                ArrowDebreu(n + 1, i) = qi * ArrowDebreu(n, i) * Df
            ElseIf n > 0 And i = n * 2 Then
                ArrowDebreu(n + 1, i) = UpProbability(n, i - 2) * ArrowDebreu(n, i - 2) * Df _
                                + (1 - UpProbability(n, i - 1) - DownProbability(n, i - 1)) * ArrowDebreu(n, i - 1) * Df _
                                + qi * ArrowDebreu(n, i) * Df
                ArrowDebreu(n + 1, i + 1) = UpProbability(n, i - 1) * ArrowDebreu(n, i - 1) * Df + (1 - pi - qi) * ArrowDebreu(n, i) * Df
                ArrowDebreu(n + 1, i + 2) = pi * ArrowDebreu(n, i) * Df
            ElseIf n > 0 And i = 1 Then
                ArrowDebreu(n + 1, i) = (1 - UpProbability(n, i - 1) - DownProbability(n, i - 1)) * ArrowDebreu(n, i - 1) * Df _
                                + qi * ArrowDebreu(n, i) * Df
            Else
            ArrowDebreu(n + 1, i) = UpProbability(n, i - 2) * ArrowDebreu(n, i - 2) * Df _
                                + (1 - UpProbability(n, i - 1) - DownProbability(n, i - 1)) * ArrowDebreu(n, i - 1) * Df _
                                + qi * ArrowDebreu(n, i) * Df
        End If
        Next
    Next
        
    If ReturnFlag = "DPM" Then
        ImpliedTrinomialTree = Application.Transpose(DownProbability)
    ElseIf ReturnFlag = "UPM" Then
        ImpliedTrinomialTree = Application.Transpose(UpProbability)
    ElseIf ReturnFlag = "DPni" Then
        ImpliedTrinomialTree = (DownProbability(STEPn, STATEi))
    ElseIf ReturnFlag = "UPni" Then
      ImpliedTrinomialTree = (UpProbability(STEPn, STATEi))
    ElseIf ReturnFlag = "ADM" Then
        ImpliedTrinomialTree = Application.Transpose(ArrowDebreu)
    ElseIf ReturnFlag = "LVM" Then
        ImpliedTrinomialTree = Application.Transpose(LocalVolatility)
    ElseIf ReturnFlag = "LVni" Then
        ImpliedTrinomialTree = Application.Transpose(LocalVolatility(STEPn, STATEi))
    ElseIf ReturnFlag = "ADni" Then
        ImpliedTrinomialTree = (ArrowDebreu(STEPn, STATEi))
    Else
    
    '// Calculation of option price using the implied trinomial tree
        If ReturnFlag = "c" Then
            z = 1
        ElseIf ReturnFlag = "p" Then
            z = -1
        End If
        For i = 0 To (2 * nSteps)
                OptionValueNode(i) = Max(0, z * (S * u ^ Max(i - nSteps, 0) * d ^ Max(nSteps - i, 0) - X))
        Next
         For n = nSteps - 1 To 0 Step -1
            For i = 0 To (n * 2)
                OptionValueNode(i) = (UpProbability(n, i) * OptionValueNode(i + 2) + (1 - UpProbability(n, i) - DownProbability(n, i)) * OptionValueNode(i + 1) + DownProbability(n, i) * OptionValueNode(i)) * Df
            Next
        Next
        ImpliedTrinomialTree = OptionValueNode(0)
    End If
End Function



'// Monte Carlo plain vanilla European option
Public Function MonteCarloStandardOption(CallPutFlag As String, S As Double, X As Double, T As Double, _
                r As Double, b As Double, v As Double, nSteps As Integer, nSimulations As Integer) As Double
            
    Dim dt As Double, St As Double
    Dim Sum As Double, Drift As Double, vSqrdt As Double
    Dim i As Integer, j As Integer, z As Integer
   
    dt = T / nSteps
    Drift = (b - v ^ 2 / 2) * dt
    vSqrdt = v * Sqr(dt)

    If CallPutFlag = "c" Then
        z = 1
    ElseIf CallPutFlag = "p" Then
        z = -1
    End If

    For i = 1 To nSimulations
        St = S
        For j = 1 To nSteps
            St = St * Exp(Drift + vSqrdt * Application.NormInv(Rnd(), 0, 1))
        Next
        Sum = Sum + Max(z * (St - X), 0)
    Next

    MonteCarloStandardOption = Exp(-r * T) * (Sum / nSimulations)

End Function

                        
'// Monte Carlo two asset Asian spread option
Public Function MonteCarloAsianSpreadOption(CallPutFlag As String, S1 As Double, S2 As Double, _
                X As Double, T As Double, r As Double, b1 As Double, b2 As Double, v1 As Double, v2 As Double, rho As Double, _
                nSteps As Integer, nSimulations As Integer) As Double
            
    Dim dt As Double, St1 As Double, St2 As Double
    Dim i As Integer, j As Integer, z As Integer
    Dim Sum As Double, Drift1 As Double, Drift2 As Double
    Dim v1Sqrdt As Double, v2Sqrdt As Double
    Dim Epsilon1 As Double, Epsilon2 As Double, Average1 As Double, Average2 As Double

    If CallPutFlag = "c" Then
        z = 1
    ElseIf CallPutFlag = "p" Then
        z = -1
    End If

    dt = T / nSteps
    Drift1 = (b1 - v1 ^ 2 / 2) * dt
    Drift2 = (b2 - v2 ^ 2 / 2) * dt
    v1Sqrdt = v1 * Sqr(dt)
    v2Sqrdt = v2 * Sqr(dt)

    For i = 1 To nSimulations
        Average1 = 0
        Average2 = 0
        St1 = S1
        St2 = S2
        For j = 1 To nSteps
            Epsilon1 = Application.NormInv(Rnd(), 0, 1)
            Epsilon2 = rho * Epsilon1 + Application.NormInv(Rnd(), 0, 1) * Sqr(1 - rho ^ 2)
            St1 = St1 * Exp(Drift1 + v1Sqrdt * Epsilon1)
            St2 = St2 * Exp(Drift2 + v2Sqrdt * Epsilon2)
            Average1 = Average1 + St1
            Average2 = Average2 + St2
        Next
        Average1 = Average1 / nSteps
        Average2 = Average2 / nSteps
        Sum = Sum + Max(z * (Average1 - Average2 - X), 0)
    Next

    MonteCarloAsianSpreadOption = Exp(-r * T) * (Sum / nSimulations)
    
End Function
