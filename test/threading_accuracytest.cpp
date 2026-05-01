/**********************************************************************
threading_accuracytest.cpp - Test that Open Babel produces identical results
when run from multiple threads simultaneously.

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
#include <mutex>
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

static mutex g_cout_mutex;

#define SAFE_COUT(msg) do { \
    lock_guard<mutex> lock(g_cout_mutex); \
    cout << msg << endl; \
} while(0)

static vector<string> loadSmiFile(const string& filename)
{
  vector<string> smiles;
  ifstream ifs(filename.c_str());
  string line;
  while (getline(ifs, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    istringstream iss(line);
    string smi;
    if (iss >> smi)
      smiles.push_back(smi);
  }
  return smiles;
}

// ----------------------------------------------------------------------
// Diagnostic: bare threads
// ----------------------------------------------------------------------
static void testBareThreads()
{
  const unsigned n = 4;
  vector<thread> threads;
  vector<int> results(n, 0);
  for (unsigned i = 0; i < n; ++i) {
    threads.emplace_back([&, i]() {
      results[i] = i * 2;
    });
  }
  for (auto& t : threads) t.join();
  for (unsigned i = 0; i < n; ++i) {
    OB_COMPARE(results[i], (int)(i * 2));
  }
  SAFE_COUT("ok 0 - bare threads");
}

// ----------------------------------------------------------------------
// Diagnostic: OBConversion per thread (no registry lookup)
// ----------------------------------------------------------------------
static void testConversionPerThread()
{
  OBConversion conv;
  OBFormat* pIn = conv.FindFormat("SMI");
  OBFormat* pOut = conv.FindFormat("CAN");
  OB_REQUIRE(pIn != nullptr);
  OB_REQUIRE(pOut != nullptr);

  const unsigned n = 4;
  vector<thread> threads;
  vector<bool> results(n, false);
  for (unsigned i = 0; i < n; ++i) {
    threads.emplace_back([&, i, pIn, pOut]() {
      OBConversion c;
      c.SetInAndOutFormats(pIn, pOut);
      OBMol mol;
      bool ok = c.ReadString(&mol, "c1ccccc1");
      results[i] = ok;
    });
  }
  for (auto& t : threads) t.join();
  for (unsigned i = 0; i < n; ++i) {
    OB_REQUIRE(results[i]);
  }
  SAFE_COUT("ok 1 - OBConversion per thread");
}

// ----------------------------------------------------------------------
// Diagnostic: AddHydrogens per thread
// ----------------------------------------------------------------------
static void testAddHydrogensPerThread()
{
  OBConversion conv;
  OBFormat* pIn = conv.FindFormat("SMI");
  OB_REQUIRE(pIn != nullptr);

  const unsigned n = 4;
  vector<thread> threads;
  vector<int> h_counts(n, 0);
  for (unsigned i = 0; i < n; ++i) {
    threads.emplace_back([&, i, pIn]() {
      OBConversion c;
      c.SetInFormat(pIn);
      OBMol mol;
      if (c.ReadString(&mol, "c1ccccc1")) {
        mol.AddHydrogens();
        h_counts[i] = mol.NumAtoms();
      }
    });
  }
  for (auto& t : threads) t.join();
  for (unsigned i = 0; i < n; ++i) {
    OB_COMPARE(h_counts[i], 12); // benzene + 6 H = 12 atoms
  }
  SAFE_COUT("ok 2 - AddHydrogens per thread");
}

// ----------------------------------------------------------------------
// Diagnostic: canonical SMILES write per thread
// ----------------------------------------------------------------------
static void testCanonicalWritePerThread()
{
  OBConversion conv;
  OBFormat* pIn = conv.FindFormat("SMI");
  OBFormat* pOut = conv.FindFormat("CAN");
  OB_REQUIRE(pIn != nullptr);
  OB_REQUIRE(pOut != nullptr);

  const unsigned n = 4;
  vector<thread> threads;
  vector<string> outputs(n);
  for (unsigned i = 0; i < n; ++i) {
    threads.emplace_back([&, i, pIn, pOut]() {
      OBConversion c;
      c.SetInAndOutFormats(pIn, pOut);
      OBMol mol;
      if (c.ReadString(&mol, "c1ccccc1")) {
        mol.AddHydrogens();
        outputs[i] = c.WriteString(&mol);
        if (!outputs[i].empty() && outputs[i].back() == '\n')
          outputs[i].pop_back();
      }
    });
  }
  for (auto& t : threads) t.join();
  for (unsigned i = 1; i < n; ++i) {
    OB_COMPARE(outputs[i], outputs[0]);
  }
  SAFE_COUT("ok 3 - canonical write per thread");
}

// ----------------------------------------------------------------------
// Part 1: SMILES roundtrip multithreaded accuracy
// ----------------------------------------------------------------------
static void testSmilesMultithreaded()
{
  string smiFile = OBTestUtil::GetFilename("aromatics.smi");
  vector<string> inputs = loadSmiFile(smiFile);
  OB_REQUIRE(!inputs.empty());
  if (inputs.size() > 200)
    inputs.resize(200);

  OBConversion conv;
  OBFormat* pIn = conv.FindFormat("SMI");
  OBFormat* pOut = conv.FindFormat("CAN");
  OB_REQUIRE(pIn != nullptr);
  OB_REQUIRE(pOut != nullptr);

  // Single-threaded baseline
  vector<string> baseline;
  for (const string& smi : inputs) {
    OBConversion c;
    c.SetInAndOutFormats(pIn, pOut);
    OBMol mol;
    OB_REQUIRE(c.ReadString(&mol, smi));
    mol.AddHydrogens();
    string can = c.WriteString(&mol);
    if (!can.empty() && can.back() == '\n')
      can.pop_back();
    baseline.push_back(can);
  }

  // Multi-threaded
  const unsigned n_threads = max(2u, thread::hardware_concurrency());
  vector<thread> threads;
  vector<vector<string>> results(n_threads);

  size_t chunk = inputs.size() / n_threads;
  size_t remainder = inputs.size() % n_threads;

  size_t start = 0;
  for (unsigned t = 0; t < n_threads; ++t) {
    size_t end = start + chunk + (t < remainder ? 1 : 0);
    threads.emplace_back([&, t, start, end, pIn, pOut]() {
      for (size_t i = start; i < end; ++i) {
        OBConversion c;
        c.SetInAndOutFormats(pIn, pOut);
        OBMol mol;
        if (!c.ReadString(&mol, inputs[i])) {
          results[t].push_back("ERROR");
          continue;
        }
        mol.AddHydrogens();
        string can = c.WriteString(&mol);
        if (!can.empty() && can.back() == '\n')
          can.pop_back();
        results[t].push_back(can);
      }
    });
    start = end;
  }

  for (auto& th : threads)
    th.join();

  // Verify
  start = 0;
  for (unsigned t = 0; t < n_threads; ++t) {
    size_t end = start + chunk + (t < remainder ? 1 : 0);
    OB_REQUIRE(results[t].size() == end - start);
    for (size_t i = start; i < end; ++i) {
      OB_COMPARE(results[t][i - start], baseline[i]);
    }
    start = end;
  }

  SAFE_COUT("ok 4 - SMILES multithreaded roundtrip");
}

// ----------------------------------------------------------------------
// Part 2: MOL2 + MMFF94 multithreaded accuracy
// ----------------------------------------------------------------------
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

static bool readMol2(const string& filename, OBMol& mol, OBFormat* pIn)
{
  ifstream ifs(filename.c_str());
  if (!ifs) return false;
  OBConversion conv;
  if (!conv.SetInFormat(pIn)) return false;
  mol.Clear();
  return conv.Read(&mol, &ifs);
}

static bool writeMol2(const OBMol& mol, string& output, OBFormat* pOut)
{
  ostringstream oss;
  OBConversion conv;
  if (!conv.SetOutFormat(pOut)) return false;
  if (!conv.Write(const_cast<OBMol*>(&mol), &oss)) return false;
  output = oss.str();
  return true;
}

static bool doubleEqual(double a, double b, double tol = 1e-4)
{
  if (std::isnan(a) && std::isnan(b)) return true;
  if (std::isnan(a) || std::isnan(b)) return false;
  return fabs(a - b) <= tol;
}

static void testMinimizeMultithreaded()
{
  string datadir = TESTDATADIR;
  vector<string> files = listMol2Files(datadir);
  if (files.empty()) {
    files.push_back(datadir + "/5sun_protein.mol2");
  }
  OB_REQUIRE(!files.empty());

  OBConversion conv;
  OBFormat* pIn = conv.FindFormat("MOL2");
  OBFormat* pOut = conv.FindFormat("MOL2");
  OB_REQUIRE(pIn != nullptr);
  OB_REQUIRE(pOut != nullptr);

  OBForceField* pFFProto = OBForceField::FindForceField("MMFF94");
  OB_REQUIRE(pFFProto != nullptr);

  // Single-threaded baseline
  vector<double> baseline_e_before, baseline_e_after;
  vector<string> baseline_mol2;
  for (const string& fn : files) {
    OBMol mol;
    OB_REQUIRE(readMol2(fn, mol, pIn));

    OBForceField* pFF = pFFProto->MakeNewInstance();
    OB_REQUIRE(pFF != nullptr);
    pFF->SetLogLevel(OBFF_LOGLVL_NONE);

    if (!pFF->Setup(mol)) {
      baseline_e_before.push_back(NAN);
      baseline_e_after.push_back(NAN);
      baseline_mol2.push_back("");
      delete pFF;
      continue;
    }

    baseline_e_before.push_back(pFF->Energy(false));
    pFF->SteepestDescent(10, 1e-6, OBFF_ANALYTICAL_GRADIENT);
    pFF->ConjugateGradients(10, 1e-6, OBFF_ANALYTICAL_GRADIENT);
    baseline_e_after.push_back(pFF->Energy(false));

    string out;
    OB_REQUIRE(writeMol2(mol, out, pOut));
    baseline_mol2.push_back(out);
    delete pFF;
  }

  // Multi-threaded
  const unsigned n_threads = max(2u, thread::hardware_concurrency());
  vector<thread> threads;
  struct ChunkResult {
    vector<double> e_before, e_after;
    vector<string> mol2;
  };
  vector<ChunkResult> results(n_threads);

  size_t chunk = files.size() / n_threads;
  size_t remainder = files.size() % n_threads;

  size_t start = 0;
  for (unsigned t = 0; t < n_threads; ++t) {
    size_t end = start + chunk + (t < remainder ? 1 : 0);
    threads.emplace_back([&, t, start, end, pIn, pOut, pFFProto]() {
      for (size_t i = start; i < end; ++i) {
        OBMol mol;
        if (!readMol2(files[i], mol, pIn)) {
          results[t].e_before.push_back(NAN);
          results[t].e_after.push_back(NAN);
          results[t].mol2.push_back("");
          continue;
        }

        OBForceField* pFF = pFFProto->MakeNewInstance();
        if (!pFF) {
          results[t].e_before.push_back(NAN);
          results[t].e_after.push_back(NAN);
          results[t].mol2.push_back("");
          continue;
        }
        pFF->SetLogLevel(OBFF_LOGLVL_NONE);

        if (!pFF->Setup(mol)) {
          results[t].e_before.push_back(NAN);
          results[t].e_after.push_back(NAN);
          results[t].mol2.push_back("");
          delete pFF;
          continue;
        }

        results[t].e_before.push_back(pFF->Energy(false));
        pFF->SteepestDescent(10, 1e-6, OBFF_ANALYTICAL_GRADIENT);
        pFF->ConjugateGradients(10, 1e-6, OBFF_ANALYTICAL_GRADIENT);
        results[t].e_after.push_back(pFF->Energy(false));

        string out;
        writeMol2(mol, out, pOut);
        results[t].mol2.push_back(out);
        delete pFF;
      }
    });
    start = end;
  }

  for (auto& th : threads)
    th.join();

  // Verify
  start = 0;
  for (unsigned t = 0; t < n_threads; ++t) {
    size_t end = start + chunk + (t < remainder ? 1 : 0);
    OB_REQUIRE(results[t].e_before.size() == end - start);
    for (size_t i = start; i < end; ++i) {
      size_t local = i - start;
      OB_ASSERT(doubleEqual(results[t].e_before[local], baseline_e_before[i]));
      OB_ASSERT(doubleEqual(results[t].e_after[local], baseline_e_after[i]));
      OB_ASSERT(results[t].mol2[local] == baseline_mol2[i]);
    }
    start = end;
  }

  SAFE_COUT("ok 5 - MMFF94 minimization multithreaded");
}

// ----------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------
int threading_accuracytest(int argc, char* argv[])
{
  int choice = 1;
  if (argc > 1)
    sscanf(argv[1], "%d", &choice);

  switch (choice) {
    case 0: testBareThreads(); break;
    case 1: testConversionPerThread(); break;
    case 2: testAddHydrogensPerThread(); break;
    case 3: testCanonicalWritePerThread(); break;
    case 4: testSmilesMultithreaded(); break;
    case 5: testMinimizeMultithreaded(); break;
    default:
      cout << "Test number " << choice << " does not exist." << endl;
      return 1;
  }
  return 0;
}
