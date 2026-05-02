/**********************************************************************
threading_performancetest.cpp - Measure multi-threaded throughput of
Open Babel's energy minimization pipeline.

Copyright (C) 2025

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
***********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <iomanip>
#include <cmath>

#include "obtest.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>

using namespace std;
using namespace OpenBabel;

#ifndef TESTDATADIR
#define TESTDATADIR="files/";
#endif

static vector<string> listMol2Files(const string& dir)
{
  vector<string> files;
  for (int i = 0; i <= 11; ++i) {
    ostringstream oss;
    oss << dir << "/culgi_" << setw(2) << setfill('0') << i << ".mol2";
    ifstream test(oss.str().c_str());
    if (test) {
      test.close();
      files.push_back(oss.str());
    }
  }
  return files;
}

static bool readMol2(const string& filename, OBMol& mol)
{
  ifstream ifs(filename.c_str());
  if (!ifs) return false;
  OBConversion conv;
  if (!conv.SetInFormat("MOL2")) return false;
  mol.Clear();
  return conv.Read(&mol, &ifs);
}

static bool writeMol2(const OBMol& mol, string& output)
{
  ostringstream oss;
  OBConversion conv;
  if (!conv.SetOutFormat("MOL2")) return false;
  if (!conv.Write(const_cast<OBMol*>(&mol), &oss)) return false;
  output = oss.str();
  return true;
}

static void runMinimizeWorkload(const vector<string>& files,
                                 int sd_steps, int cg_steps,
                                 vector<double>& out_energies)
{
  OBForceField* pFFProto = OBForceField::FindForceField("MMFF94");
  if (!pFFProto) {
    cerr << "MMFF94 not found" << endl;
    return;
  }

  for (const string& fn : files) {
    OBMol mol;
    if (!readMol2(fn, mol)) continue;

    OBForceField* pFF = pFFProto->MakeNewInstance();
    if (!pFF) continue;
    pFF->SetLogLevel(OBFF_LOGLVL_NONE);

    if (!pFF->Setup(mol)) {
      delete pFF;
      continue;
    }

    pFF->Energy(false);
    if (sd_steps > 0)
      pFF->SteepestDescent(sd_steps, 1e-6, OBFF_ANALYTICAL_GRADIENT);
    if (cg_steps > 0)
      pFF->ConjugateGradients(cg_steps, 1e-6, OBFF_ANALYTICAL_GRADIENT);
    double e = pFF->Energy(false);
    out_energies.push_back(e);

    string out;
    writeMol2(mol, out);

    delete pFF;
  }
}

static double benchmark(const vector<string>& files,
                        unsigned n_threads,
                        int sd_steps, int cg_steps,
                        int repeat)
{
  // Warmup
  {
    vector<double> energies;
    runMinimizeWorkload(files, sd_steps, cg_steps, energies);
  }

  auto start = chrono::steady_clock::now();

  for (int r = 0; r < repeat; ++r) {
    if (n_threads <= 1) {
      vector<double> energies;
      runMinimizeWorkload(files, sd_steps, cg_steps, energies);
    } else {
      vector<thread> threads;
      vector<vector<double>> thread_energies(n_threads);

      size_t chunk = files.size() / n_threads;
      size_t remainder = files.size() % n_threads;
      size_t start_idx = 0;

      for (unsigned t = 0; t < n_threads; ++t) {
        size_t end_idx = start_idx + chunk + (t < remainder ? 1 : 0);
        vector<string> chunkFiles(files.begin() + start_idx, files.begin() + end_idx);
        threads.emplace_back([&, t, chunkFiles]() mutable {
          vector<double> e;
          runMinimizeWorkload(chunkFiles, sd_steps, cg_steps, e);
          thread_energies[t] = std::move(e);
        });
        start_idx = end_idx;
      }

      for (auto& th : threads)
        th.join();
    }
  }

  auto end = chrono::steady_clock::now();
  chrono::duration<double> elapsed = end - start;
  return elapsed.count();
}

int threading_performancetest(int argc, char* argv[])
{
  string datadir = TESTDATADIR;
  vector<string> files = listMol2Files(datadir);
  if (files.empty()) {
    cerr << "No MOL2 test files found" << endl;
    return 1;
  }

  int choice = 1;
  if (argc > 1)
    sscanf(argv[1], "%d", &choice);

  int sd_steps = 10;
  int cg_steps = 10;
  int repeat = 10; // repeat the whole set to get measurable times

  if (choice == 2) {
    // Quick sanity check mode (1 repeat)
    repeat = 1;
  }

  cout << "# Open Babel Threading Performance Test" << endl;
  cout << "# Molecules: " << files.size() << endl;
  cout << "# SD steps: " << sd_steps << ", CG steps: " << cg_steps << endl;
  cout << "# Repeat: " << repeat << endl;
  cout << "#" << endl;
  cout << "# threads  time(s)  molecules  mol/sec  speedup  efficiency" << endl;

  double baseline_time = 0.0;
  unsigned max_threads = max(1u, thread::hardware_concurrency());
  vector<unsigned> thread_counts = {1};
  for (unsigned t = 2; t <= max_threads; t *= 2)
    thread_counts.push_back(t);
  if (thread_counts.back() != max_threads && max_threads > 1)
    thread_counts.push_back(max_threads);

  for (unsigned n_threads : thread_counts) {
    double t = benchmark(files, n_threads, sd_steps, cg_steps, repeat);
    size_t total_mols = files.size() * repeat;
    double throughput = total_mols / t;

    if (n_threads == 1)
      baseline_time = t;

    double speedup = baseline_time > 0 ? (baseline_time / t) : 1.0;
    double efficiency = speedup / n_threads;

    cout << setw(9) << n_threads
         << "  " << fixed << setprecision(3) << setw(7) << t
         << "  " << setw(9) << total_mols
         << "  " << fixed << setprecision(1) << setw(7) << throughput
         << "  " << fixed << setprecision(2) << speedup
         << "  " << fixed << setprecision(2) << efficiency
         << endl;
  }

  return 0;
}
