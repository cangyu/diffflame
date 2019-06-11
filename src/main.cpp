#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/Inlet1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

using namespace Cantera;
using namespace std;

const string DATA_DIR = "../data/";
const bool USE_EXISTING_DATA = true;
const size_t INITIAL_PNT_NUM = 11;
const size_t DATA_WIDTH = 32;
const size_t DATA_DIGITS = 18;

template<typename T>
T relaxation(const T& a, const T& b, double alpha)
{
    return (1-alpha) * a + alpha * b;
}

double calc_density(double P, double T, const double *Y, const double *MW, size_t K)
{
    double ret = 0.0;
    for(auto k = 0; k < K; ++k)
        ret += Y[k] / MW[k];
    ret = P/(GasConstant * T * ret);
    return ret;
}

void diffflame(double mdot_f, double mdot_o, double domain_length)
{
    // Use the GRI3.0 mechanism
    IdealGasMix gas("gri30.cti", "gri30");
    const auto K = gas.nSpecies();
    cout << K << " species" << endl;
    const auto MW = gas.molecularWeights();
    const auto NAME = gas.speciesNames();

    // Constant pressure and input mass flux
    const double P = OneAtm;
    const double T_o = 300.0;
    const double T_f = 300.0;
    
    // Initial grid
    const size_t N = INITIAL_PNT_NUM;
    const double dz = domain_length/(N-1);
    vector<double> z(N, 0.0);
    z[0] = -domain_length / 2;
    for(auto i = 1; i < N; i++)
        z[i] = z[i-1] + dz;

    // Inlet boundaries of flow field
    Inlet1D fuel_inlet;
    fuel_inlet.setMoleFractions("CH4:0.5, H2:0.5");
    fuel_inlet.setMdot(mdot_f);
    fuel_inlet.setTemperature(T_f);

    Inlet1D oxidizer_inlet;
    oxidizer_inlet.setMoleFractions("N2:0.78, O2:0.21, AR:0.01");
    oxidizer_inlet.setMdot(mdot_o);
    oxidizer_inlet.setTemperature(T_o);

    // Interior of flow field
    const vector<double> tol_ss{1.0e-5, 1.0e-9}; // [rtol atol] for steady-state
    const vector<double> tol_ts{1.0e-3, 1.0e-9}; // [rtol atol] for time-stepping
    StFlow flow(&gas, K, N);
    flow.setAxisymmetricFlow();
    flow.setupGrid(N, z.data());
    Transport* tr = newTransportMgr("UnityLewis", &gas);
    flow.setTransport(*tr);
    flow.setKinetics(gas);
    flow.setSteadyTolerances(tol_ss[0], tol_ss[1]);
    flow.setTransientTolerances(tol_ts[0], tol_ts[1]);

    // The simulation domain(for calculating the residuals)
    vector<Domain1D*> domains{&fuel_inlet, &flow, &oxidizer_inlet};
    Sim1D flame(domains);
    
    // Init
    const auto Teq = relaxation(T_f, T_o, 0.5);
    gas.setState_TPX(Teq, P, "CH4:1.0, H2:1.0, N2:9.2857, O2: 2.5, AR:0.1190");
    try {
        gas.equilibrate("HP");
    } catch (CanteraError& err) {
        cout << err.what() << endl;
    }
    const auto Tad = gas.temperature();
    cout << "Tad = " << Tad << "K" << endl;
    vector<double> Xst(K), Yst(K);
    gas.getMoleFractions(Xst.data());
    gas.getMassFractions(Yst.data());

    vector<double> Y_f(K), Y_o(K);
    for (auto k = 0; k < K; ++k)
    {
        Y_f[k] = fuel_inlet.massFraction(k);
        Y_o[k] = oxidizer_inlet.massFraction(k);
    }
    const double rho_f = calc_density(P, T_f, Y_f.data(), MW.data(), K);
    const double u_f = mdot_f / rho_f;
    const double rho_o = calc_density(P, T_o, Y_o.data(), MW.data(), K);
    const double u_o = -mdot_o / rho_o;

    vector_fp locs{0, 0.0, 1.0};
    vector_fp values{T_f, Tad, T_o};
    flame.setInitialGuess("T", locs, values);

    for(auto k = 0; k < K; ++k)
    {
        values[0] = Y_f[k];
        values[1] = Yst[k];
        values[2] = Y_o[k];
        flame.setInitialGuess(NAME[k],locs,values);
    }

    // Report B.C.
    cout << "Domain length: " <<  domain_length << " m" << endl;
    cout << "B.C. at Fuel inlet:" << endl;
    cout << "    Density: " << rho_f << " Kg/m^3" << endl;
    cout << "    Velocity: " << u_f << " m/s" << endl;
    fuel_inlet.showSolution(nullptr);
    cout << "B.C. at Oxidizer inlet:" << endl;
    cout << "    Density: " << rho_o << " Kg/m^3" << endl;
    cout << "    Velocity: " << u_o << " m/s" << endl;
    oxidizer_inlet.showSolution(nullptr);

    // Solve
    const int IdxOfFlowDomain = 1;
    const int MaxNumOfPnt = 501;
    flame.setMaxGridPoints(IdxOfFlowDomain,MaxNumOfPnt);
    const double ratio=200.0;
    const double slope=0.1;
    const double curve=0.2;
    const double prune=-0.1;
    flame.setRefineCriteria(IdxOfFlowDomain,ratio,slope,curve,prune);   
    flow.setPressure(P);
    flow.solveEnergyEqn(); 
    flame.solve(1,true);

    // Output
    stringstream ss;
    ss << "mf=" << mdot_f << "_mo=" << mdot_o << "_L=" << domain_length << "_raw.txt";
    ofstream fout(ss.str());
    
    fout << setw(DATA_WIDTH) << "x" << setw(DATA_WIDTH) << "u" << setw(DATA_WIDTH) << "V" << setw(DATA_WIDTH) << "T" << setw(DATA_WIDTH) << "Lambda";
    for(auto k = 0; k < K; ++k)
    {
        string yk_name = "Y_" + NAME[k];
        fout << setw(DATA_WIDTH) << yk_name;
    }
    fout << endl;

    double Tmax = 0.0;
    for(auto i = 0; i < flow.nPoints(); ++i)
    {
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flow.grid(i);
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 0, i); // u
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 1, i); // V
        double loc_T = flame.value(IdxOfFlowDomain, 2, i);
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << loc_T; // T
        if(loc_T > Tmax)
            Tmax = loc_T;
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 3, i); // Lambda
        for(auto k = 0; k < K; ++k)
        {
            auto yk = flame.value(IdxOfFlowDomain, 5+k, i);
            if(yk < 0.0)
                yk = 0.0;
            if(yk > 1.0)
                yk = 1.0;
            fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << yk; // Yk
        }
        fout << endl;
    }
    fout.close();

    // Report 
    cout << "Tmax=" << Tmax << "K" << endl;
}

void diffflame(double mdot_f, double mdot_o, double domain_length, const vector<vector<double>> &init_data)
{
    // Setup initial grid
    const size_t N = init_data.size();
    cout << "Initialize using existing data with " << N << " points!" << endl;
    vector<double> z(N, 0.0);
    for(auto i = 0; i < N; i++)
        z[i] = init_data[i][0];

    // Use the GRI3.0 mechanism
    IdealGasMix gas("gri30.cti", "gri30");
    const auto K = gas.nSpecies();
    cout << K << " species" << endl;
    const auto MW = gas.molecularWeights();
    const auto NAME = gas.speciesNames();

    // Constant pressure and input mass flux
    const double P = OneAtm;
    const double T_o = 300.0;
    const double T_f = 300.0;

    // Inlet boundaries of flow field
    Inlet1D fuel_inlet;
    fuel_inlet.setMoleFractions("CH4:0.5, H2:0.5");
    fuel_inlet.setMdot(mdot_f);
    fuel_inlet.setTemperature(T_f);

    Inlet1D oxidizer_inlet;
    oxidizer_inlet.setMoleFractions("N2:0.78, O2:0.21, AR:0.01");
    oxidizer_inlet.setMdot(mdot_o);
    oxidizer_inlet.setTemperature(T_o);

    // Interior of flow field
    const vector<double> tol_ss{1.0e-5, 1.0e-9}; // [rtol atol] for steady-state
    const vector<double> tol_ts{1.0e-3, 1.0e-9}; // [rtol atol] for time-stepping
    StFlow flow(&gas, K, N);
    flow.setAxisymmetricFlow();
    flow.setupGrid(N, z.data());
    Transport* tr = newTransportMgr("UnityLewis", &gas);
    flow.setTransport(*tr);
    flow.setKinetics(gas);
    flow.setSteadyTolerances(tol_ss[0], tol_ss[1]);
    flow.setTransientTolerances(tol_ts[0], tol_ts[1]);

    // The simulation domain(for calculating the residuals)
    vector<Domain1D*> domains{&fuel_inlet, &flow, &oxidizer_inlet};
    Sim1D flame(domains);
    
    // Init
    const auto Teq = relaxation(T_f, T_o, 0.5);
    gas.setState_TPX(Teq, P, "CH4:1.0, H2:1.0, N2:9.2857, O2: 2.5, AR:0.1190");
    try {
        gas.equilibrate("HP");
    } catch (CanteraError& err) {
        cout << err.what() << endl;
    }
    const auto Tad = gas.temperature();
    cout << "Tad = " << Tad << "K" << endl;
    vector<double> Xst(K), Yst(K);
    gas.getMoleFractions(Xst.data());
    gas.getMassFractions(Yst.data());

    vector<double> Y_f(K), Y_o(K);
    for (auto k = 0; k < K; ++k)
    {
        Y_f[k] = fuel_inlet.massFraction(k);
        Y_o[k] = oxidizer_inlet.massFraction(k);
    }
    const double rho_f = calc_density(P, T_f, Y_f.data(), MW.data(), K);
    const double u_f = mdot_f / rho_f;
    const double rho_o = calc_density(P, T_o, Y_o.data(), MW.data(), K);
    const double u_o = -mdot_o / rho_o;

    vector_fp locs;
    for(int i = 0; i < init_data.size(); ++i)
    {
        double cur_pos = (init_data[i][0] - init_data[0][0])/domain_length;
        locs.push_back(cur_pos);
    }

    vector_fp values(locs.size(), 0.0);
    for(int i = 0; i < init_data.size(); ++i)
        values[i] = init_data[i][1];
    flame.setInitialGuess("u", locs, values);
    
    for(int i = 0; i < init_data.size(); ++i)
        values[i] = init_data[i][2];
    flame.setInitialGuess("V", locs, values);

    for(int i = 0; i < init_data.size(); ++i)
        values[i] = init_data[i][3];
    flame.setInitialGuess("T", locs, values);

    for(int i = 0; i < init_data.size(); ++i)
        values[i] = init_data[i][4];
    flame.setInitialGuess("lambda", locs, values);
    
    for(auto k = 0; k < K; ++k)
    {
        for(int i = 0; i < init_data.size(); ++i)
            values[i] = init_data[i][5+k];
        flame.setInitialGuess(NAME[k],locs,values);
    }

    // Report B.C.
    cout << "Domain length: " <<  domain_length << " m" << endl;
    cout << "B.C. at Fuel inlet:" << endl;
    cout << "    Density: " << rho_f << " Kg/m^3" << endl;
    cout << "    Velocity: " << u_f << " m/s" << endl;
    fuel_inlet.showSolution(nullptr);
    cout << "B.C. at Oxidizer inlet:" << endl;
    cout << "    Density: " << rho_o << " Kg/m^3" << endl;
    cout << "    Velocity: " << u_o << " m/s" << endl;
    oxidizer_inlet.showSolution(nullptr);

    // Solve
    const int IdxOfFlowDomain = 1;
    const int MaxNumOfPnt = 501;
    flame.setMaxGridPoints(IdxOfFlowDomain,MaxNumOfPnt);
    const double ratio=200.0;
    const double slope=0.1;
    const double curve=0.2;
    const double prune=-0.1;
    flame.setRefineCriteria(IdxOfFlowDomain,ratio,slope,curve,prune);   
    flow.setPressure(P);
    flow.solveEnergyEqn(); 
    flame.solve(1,true);

    // Output
    stringstream ss;
    ss << "mf=" << mdot_f << "_mo=" << mdot_o << "_L=" << domain_length << "_raw.txt";
    ofstream fout(ss.str());
    
    fout << setw(DATA_WIDTH) << "x" << setw(DATA_WIDTH) << "u" << setw(DATA_WIDTH) << "V" << setw(DATA_WIDTH) << "T" << setw(DATA_WIDTH) << "Lambda";
    for(auto k = 0; k < K; ++k)
    {
        string yk_name = "Y_" + NAME[k];
        fout << setw(DATA_WIDTH) << yk_name;
    }
    fout << endl;

    double Tmax = 0.0;
    for(auto i = 0; i < flow.nPoints(); ++i)
    {
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flow.grid(i);
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 0, i); // u
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 1, i); // V
        double loc_T = flame.value(IdxOfFlowDomain, 2, i);
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << loc_T; // T
        if(loc_T > Tmax)
            Tmax = loc_T;
        fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << flame.value(IdxOfFlowDomain, 3, i); // Lambda
        for(auto k = 0; k < K; ++k)
        {
            auto yk = flame.value(IdxOfFlowDomain, 5+k, i);
            if(yk < 0.0)
                yk = 0.0;
            if(yk > 1.0)
                yk = 1.0;
            fout << setw(DATA_WIDTH) << setprecision(DATA_DIGITS) << scientific << yk; // Yk
        }
        fout << endl;
    }
    fout.close();

    // Report 
    cout << "Tmax=" << Tmax << "K" << endl;
}

int main(int argc, char *argv[])
{   
    double mf=0.0, mo=0.0, L=0.1;
    
    // Get input parameters
    if(argc == 4)
    {
        mf = atof(argv[1]);
        mo = atof(argv[2]);
        L = atof(argv[3]);
    }
    else
        return -1;

    // Load existing data
    stringstream ss;
    ss << DATA_DIR << "mf=" << mf << "_mo=" << mo << "_L=" << L << "_raw.txt";
    ifstream fin(ss.str());
    bool init_data_available = fin.good();
    fin.close();
    vector<vector<double>> data;
    if(init_data_available)
    {
        char c = 0;
        int NumOfLine = 0;

        fin.open(ss.str());
        while(fin.get(c))
        {
            if (c=='\n')
                ++NumOfLine;
        }
        assert(fin.eof());
        fin.close();
        
        fin.open(ss.str());
        while(fin.get(c))
        {
            if(c=='\n')
                break;
        }

        vector<double> cur_data(58, 0.0);
        for(int line = 0; line<NumOfLine-1; ++line)
        {
            for(int k = 0; k < 58; ++k)
                fin>>cur_data[k];
            data.push_back(cur_data);
        }
        fin.close();
    }
    
    // Solve
    try {
        if(init_data_available && USE_EXISTING_DATA)
            diffflame(mf, mo, L, data);
        else
            diffflame(mf, mo, L);
    }
    catch (CanteraError& err) {
        cout << err.what() << endl;
    }
    return 0;
}
