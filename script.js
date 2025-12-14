/* ==========================================================================
   TUGAS BESAR OS2103: KOMPUTASI OSEANOGRAFI (WEB VERSION)
   Implementasi Lengkap Algoritma UNESCO No. 44 (1983)
   ========================================================================== */

// --- KONSTANTA GLOBAL ---
const k_sal = 0.0162;

// ==========================================================================
// BAGIAN 1: RUMUS MATEMATIKA (ALGORITMA 1-9)
// ==========================================================================

// 1. CONDUCTIVITY -> SALINITY (PSS-78)
function algo1_salinity(R, T, P) {
    const a0=0.0080, a1=-0.1692, a2=25.3851, a3=14.0941, a4=-7.0261, a5=2.7081;
    const b0=0.0005, b1=-0.0056, b2=-0.0066, b3=-0.0375, b4=0.0636, b5=-0.0144;
    const c0=0.6766097, c1=2.00564E-2, c2=1.104259E-4, c3=-6.9698E-7, c4=1.0031E-9;
    const d1=3.426E-2, d2=4.464E-4, d3=4.215E-1, d4=-3.107E-3;
    const e1=2.070E-5, e2=-6.370E-10, e3=3.989E-15;

    let DT = T - 15.0;
    let Rp = 1.0 + (P * (e1 + e2*P + e3*Math.pow(P,2))) / (1.0 + d1*T + d2*Math.pow(T,2) + (d3 + d4*T)*R);
    let rt = c0 + T*(c1 + T*(c2 + T*(c3 + T*c4)));
    let Rt = R / (Rp * rt);
    let RT = Math.sqrt(Math.abs(Rt));

    let S = a0 + RT*(a1 + RT*(a2 + RT*(a3 + RT*(a4 + RT*a5)))) + 
            (DT / (1.0 + k_sal*DT)) * (b0 + RT*(b1 + RT*(b2 + RT*(b3 + RT*(b4 + RT*b5)))));
    return Math.max(0, S);
}

// 2. SALINITY -> CONDUCTIVITY (Iterative Newton-Raphson)
function algo2_conductivity(S, T, P) {
    const a1=-0.1692, a2=25.3851, a3=14.0941, a4=-7.0261, a5=2.7081;
    const b1=-0.0056, b2=-0.0066, b3=-0.0375, b4=0.0636, b5=-0.0144;
    const c0=0.6766097, c1=2.00564E-2, c2=1.104259E-4, c3=-6.9698E-7, c4=1.0031E-9;
    const d1=3.426E-2, d2=4.464E-4, d3=4.215E-1, d4=-3.107E-3;
    const e1=2.070E-5, e2=-6.370E-10, e3=3.989E-15;

    let DT = T - 15.0;
    let RT = Math.sqrt(S / 35.0);
    let DELS = 1.0;
    let i = 0;

    // Iterasi Newton-Raphson
    while (DELS > 1.0E-5 && i < 20) {
        i++;
        let SI = 0.0080 + RT*(a1 + RT*(a2 + RT*(a3 + RT*(a4 + RT*a5)))) + 
                 (DT/(1.0+k_sal*DT)) * (0.0005 + RT*(b1 + RT*(b2 + RT*(b3 + RT*(b4 + RT*b5)))));
        let DERIV = a1 + RT*(2*a2 + RT*(3*a3 + RT*(4*a4 + RT*5*a5))) + 
                    (DT/(1.0+k_sal*DT)) * (b1 + RT*(2*b2 + RT*(3*b3 + RT*(4*b4 + RT*5*b5))));
        RT = RT + (S - SI) / DERIV;
        DELS = Math.abs(S - SI);
    }

    let rt = c0 + T*(c1 + T*(c2 + T*(c3 + T*c4)));
    let A = d3 + d4*T;
    let B = 1.0 + d1*T + d2*Math.pow(T,2);
    let C = P * (e1 + e2*P + e3*Math.pow(P,2));
    let R_pure = Math.pow(RT, 2) * rt;
    let Rp = 1.0 + C / (B + A * R_pure);
    
    return R_pure * Rp;
}

// 3. DENSITY (EOS-80)
function algo3_density(S, T, P) {
    let P_bar = P / 10.0;
    let rho_w = 999.842594 + 6.793952E-2*T - 9.095290E-3*Math.pow(T,2) + 
                1.001685E-4*Math.pow(T,3) - 1.120083E-6*Math.pow(T,4) + 6.536332E-9*Math.pow(T,5);
    
    let K0 = rho_w + (0.824493 - 0.0040899*T + 7.6438E-5*Math.pow(T,2) - 8.2467E-7*Math.pow(T,3) + 5.3875E-9*Math.pow(T,4))*S + 
             (-0.00572466 + 1.0227E-4*T - 1.6546E-6*Math.pow(T,2))*S*Math.sqrt(S) + 4.8314E-4*Math.pow(S,2);

    let Kw = 19652.21 + 148.4206*T - 2.327105*Math.pow(T,2) + 1.360477E-2*Math.pow(T,3) - 5.155288E-5*Math.pow(T,4);
    let BulkMod = Kw + (54.6746 - 0.603459*T + 1.09987E-2*Math.pow(T,2) - 6.1670E-5*Math.pow(T,3))*S + 
                  (0.07944 + 0.016483*T - 5.3009E-4*Math.pow(T,2))*S*Math.sqrt(S);
    
    let A = 3.239908 + 1.43713E-3*T + 1.16092E-4*Math.pow(T,2) - 5.77905E-7*Math.pow(T,3) + 
            (2.2838E-3 - 1.0981E-5*T - 1.6078E-6*Math.pow(T,2))*S + 1.91075E-4*S*Math.sqrt(S);
    let B = 8.50935E-5 - 6.12293E-6*T + 5.2787E-8*Math.pow(T,2) + (-9.9348E-7 + 2.0816E-8*T + 9.1697E-10*Math.pow(T,2))*S;

    BulkMod = BulkMod + A*P_bar + B*Math.pow(P_bar,2);
    return K0 / (1.0 - P_bar / BulkMod);
}

// 4. PRESSURE TO DEPTH
function algo4_depth(P, LAT) {
    let rad = LAT / 57.29578; // Deg to Rad
    let X = Math.pow(Math.sin(rad), 2);
    let GR = 9.780318 * (1.0 + (5.2788E-3 + 2.36E-5*X)*X) + 1.092E-6*P;
    let Z = (((-1.82E-15*P + 2.279E-10)*P - 2.2512E-5)*P + 9.72659)*P;
    return Z / GR;
}

// 5. FREEZING POINT
function algo5_freezing(S, P) {
    const a0 = -0.0575, a1 = 1.710523E-3, a2 = -2.154996E-4, b = -7.53E-4;
    return (a0*S) + (a1*S*Math.sqrt(S)) + (a2*Math.pow(S,2)) + (b*P);
}

// 6. SPECIFIC HEAT
function algo6_specific_heat(S, T, P) {
    let P_bar = P / 10.0;
    let SR = Math.sqrt(S);
    
    // Cp at Atmospheric Pressure
    let Cp_0 = 4217.4 - 3.720283*T + 0.1412855*Math.pow(T,2) - 2.654387E-3*Math.pow(T,3) + 2.093236E-5*Math.pow(T,4) +
               (-7.643575 + 0.1072763*T - 1.38385E-3*Math.pow(T,2))*S +
               (0.1770383 - 4.07718E-3*T + 5.148E-5*Math.pow(T,2))*S*SR;
    
    // Pressure terms
    let A = -0.49592 + 1.45747E-2*T - 3.13885E-4*Math.pow(T,2) + 2.0357E-6*Math.pow(T,3) + 1.7168E-8*Math.pow(T,4);
    let B = 2.4931E-4 - 1.08645E-5*T + 2.87533E-7*Math.pow(T,2) - 4.0027E-9*Math.pow(T,3) + 2.2956E-11*Math.pow(T,4);
    let C = -5.422E-8 + 2.6380E-9*T - 6.5637E-11*Math.pow(T,2) + 6.136E-13*Math.pow(T,3);
    let Del_Cp1 = A*P_bar + B*Math.pow(P_bar,2) + C*Math.pow(P_bar,3);

    A = (4.9247E-3 - 1.28315E-4*T + 9.802E-7*Math.pow(T,2) + 2.5941E-8*Math.pow(T,3) - 2.9179E-10*Math.pow(T,4)) * S +
        (-1.2331E-4 - 1.517E-6*T + 3.122E-8*Math.pow(T,2)) * S*SR;
    B = (-2.9558E-6 + 1.17054E-7*T - 2.3905E-9*Math.pow(T,2) + 1.8448E-11*Math.pow(T,3)) * S +
        (9.971E-8) * S*SR;
    C = (5.540E-10 - 1.7682E-11*T + 3.513E-13*Math.pow(T,2)) * S + (-1.4300E-12 * T) * S*SR;
    let Del_Cp2 = A*P_bar + B*Math.pow(P_bar,2) + C*Math.pow(P_bar,3);

    return Cp_0 + Del_Cp1 + Del_Cp2;
}

// 7. ADIABATIC LAPSE RATE
function algo7_adiabatic(S, T, P) {
    let DS = S - 35.0;
    return 3.5803E-5 + 8.5258E-6*T - 6.8360E-8*Math.pow(T,2) + 6.6228E-10*Math.pow(T,3) +
           (1.8932E-6 - 4.2393E-8*T)*DS +
           (1.8741E-8 - 6.7795E-10*T + 8.7330E-12*Math.pow(T,2) - 5.4481E-14*Math.pow(T,3))*P +
           (-1.1351E-10 + 2.7759E-12*T)*DS*P +
           (-4.6206E-13 + 1.8676E-14*T - 2.1687E-16*Math.pow(T,2))*Math.pow(P,2);
}

// 8. POTENTIAL TEMPERATURE (Runge-Kutta 4th Order Integration)
function algo8_potential_temp(S, T, P, P_ref = 0) {
    let H = P_ref - P;
    let T_local = T;
    
    // RK4 Step 1
    let Q1 = algo7_adiabatic(S, T_local, P);
    // RK4 Step 2
    let Q2 = algo7_adiabatic(S, T_local + 0.5*H*Q1, P + 0.5*H);
    // RK4 Step 3
    let Q3 = algo7_adiabatic(S, T_local + 0.5*H*Q2, P + 0.5*H);
    // RK4 Step 4
    let Q4 = algo7_adiabatic(S, T_local + H*Q3, P + H);

    return T_local + (H/6.0) * (Q1 + 2.0*Q2 + 2.0*Q3 + Q4);
}

// 9. SOUND SPEED
function algo9_sound(S, T, P) {
    let Pb = P / 10.0; // Bar
    let SR = Math.sqrt(Math.abs(S));
    let Cw = 1402.388 + 5.03711*T - 0.0580852*Math.pow(T,2) + 3.3420E-4*Math.pow(T,3) - 
             1.47800E-6*Math.pow(T,4) + 3.1464E-9*Math.pow(T,5) +
             (0.153563 + 6.8982E-4*T - 8.1788E-6*Math.pow(T,2) + 1.3621E-7*Math.pow(T,3) - 6.1185E-10*Math.pow(T,4))*Pb +
             (3.1260E-5 - 1.7107E-6*T + 2.5974E-8*Math.pow(T,2) - 2.5335E-10*Math.pow(T,3) + 1.0405E-12*Math.pow(T,4))*Math.pow(Pb,2) +
             (-9.7729E-9 + 3.8504E-10*T - 2.3643E-12*Math.pow(T,2))*Math.pow(Pb,3);

    let A = 1.389 - 1.262E-2*T + 7.164E-5*Math.pow(T,2) + 2.006E-6*Math.pow(T,3) - 3.21E-8*Math.pow(T,4) +
            (9.4742E-5 - 1.2580E-5*T - 6.4885E-8*Math.pow(T,2) + 1.0507E-8*Math.pow(T,3) - 2.0122E-10*Math.pow(T,4))*Pb +
            (-3.9064E-7 + 9.1041E-9*T - 1.6002E-10*Math.pow(T,2) + 7.988E-12*Math.pow(T,3))*Math.pow(Pb,2) +
            (1.100E-10 + 6.649E-12*T - 3.389E-13*Math.pow(T,2))*Math.pow(Pb,3);

    let B = -1.922E-2 - 4.42E-5*T + (7.3637E-5 + 1.7945E-7*T)*Pb;
    let D = 1.727E-3 - 7.9836E-6*Pb;
    return Cw + A*S + B*S*SR + D*Math.pow(S,2);
}

// ==========================================================================
// BAGIAN 2: LOGIKA UI (User Interface)
// ==========================================================================

function updateInputs() {
    let choice = document.getElementById("algoSelect").value;
    
    // Ambil elemen Label dan Input
    let l1 = document.getElementById("val1").previousElementSibling;
    let l2 = document.getElementById("val2").previousElementSibling;
    let l3 = document.getElementById("val3").previousElementSibling; // Label Pressure
    
    let i1 = document.getElementById("val1");
    let i2 = document.getElementById("val2");
    let i3 = document.getElementById("val3");

    // --- RESET DEFAULT (Untuk Algo 1, 2, 3, 6, 7, 9) ---
    l1.textContent = "Salinity (PSU) / Ratio:";
    l2.textContent = "Temperature (C):";
    l3.textContent = "Pressure (decibars):";
    i1.style.display = "block"; l1.style.display = "block";
    i2.style.display = "block"; l2.style.display = "block";
    i3.style.display = "block"; l3.style.display = "block";

    // --- KUSTOMISASI PER ALGORITMA ---
    if (choice == "4") { // Pressure to Depth
        l1.textContent = "Pressure (decibars):";
        l2.textContent = "Latitude (degrees):";
        i3.style.display = "none"; l3.style.display = "none"; // Hide input ke-3
    } 
    else if (choice == "5") { // Freezing Point
        l1.textContent = "Salinity (PSU):";
        // Hide Temp
        i2.style.display = "none"; l2.style.display = "none"; 
        // Input ke-3 jadi Pressure
        l3.textContent = "Pressure (decibars):";
    }
    else if (choice == "8") { // Potential Temp
        // Butuh P_Ref (Reference Pressure). Kita pakai slot 3 input aja biar simpel.
        // Asumsi P_Ref = 0 (Permukaan) kalau tidak ada input ke-4.
        l3.textContent = "Pressure In-Situ (db):";
    }
}

function calculate() {
    let choice = document.getElementById("algoSelect").value;
    let v1 = parseFloat(document.getElementById("val1").value);
    let v2 = parseFloat(document.getElementById("val2").value);
    let v3 = parseFloat(document.getElementById("val3").value);
    
    if (isNaN(v1)) v1 = 0; // Safety check
    if (isNaN(v2)) v2 = 0;
    if (isNaN(v3)) v3 = 0;

    let hasil = 0;
    let unit = "";

    switch(choice) {
        case "1":
            hasil = algo1_salinity(v1, v2, v3);
            unit = " PSU";
            break;
        case "2":
            hasil = algo2_conductivity(v1, v2, v3);
            unit = " (Ratio)";
            break;
        case "3":
            hasil = algo3_density(v1, v2, v3);
            unit = " kg/m³";
            break;
        case "4": // v1=Press, v2=Lat
            hasil = algo4_depth(v1, v2);
            unit = " meter";
            break;
        case "5": // v1=Sal, v3=Press (v2 hidden)
            hasil = algo5_freezing(v1, v3);
            unit = " °C";
            break;
        case "6":
            hasil = algo6_specific_heat(v1, v2, v3);
            unit = " J/(kg °C)";
            break;
        case "7":
            hasil = algo7_adiabatic(v1, v2, v3);
            unit = " °C/dbar";
            break;
        case "8": // Pot Temp (Ref=0)
            hasil = algo8_potential_temp(v1, v2, v3, 0.0);
            unit = " °C (at Surface)";
            break;
        case "9":
            hasil = algo9_sound(v1, v2, v3);
            unit = " m/s";
            break;
    }

    // Tampilkan Hasil
    document.getElementById("result").style.display = "block";
    document.getElementById("outputValue").innerHTML = hasil.toFixed(5) + unit;
}

// Jalankan updateInputs() sekali saat halaman dimuat
window.onload = updateInputs;
