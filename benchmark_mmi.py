#!/usr/bin/env python3
"""Benchmark MMFF94 energy calculation overhead of MakeNewInstance vs prototype reuse."""

import time
import sys
sys.path.insert(0, '/Users/jt/PycharmProjects/openbabel/scripts/python')

from openbabel import openbabel as ob
from openbabel import pybel

# Test molecules
smiles_list = [
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
    "c1ccc2c(c1)c(c[nH]2)CCN",      # tryptamine
    "CC(=O)Oc1ccccc1C(=O)O",        # aspirin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # caffeine
    "c1cc(ccc1C(=O)O)O",            # salicylic acid
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "C1=CC=C(C=C1)CC(C(=O)O)N",
    "C(C(=O)O)N",
    "CC(C)C1=CC=C(C=C1)C(C)C(=O)O",
    "C1=CC=C2C(=C1)C(=CN2)CCN",
]

# Warm up - find forcefield
ff_proto = ob.OBForceField.FindForceField("MMFF94")
if not ff_proto:
    print("MMFF94 not found")
    sys.exit(1)

# Parse molecules
mols = []
for smi in smiles_list:
    mol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.ReadString(mol, smi)
    mols.append(mol)

N = 1000

# Benchmark 1: MakeNewInstance (our current code path)
print("Benchmarking MakeNewInstance + Setup + Energy...")
t0 = time.time()
for i in range(N):
    mol = mols[i % len(mols)]
    ff = ff_proto.MakeNewInstance()
    ff.Setup(mol)
    ff.Energy()
    del ff
t1 = time.time()
rate_newinstance = N / (t1 - t0)
print(f"  MakeNewInstance: {N} runs in {t1-t0:.3f}s = {rate_newinstance:.1f} mol/s")

# Benchmark 2: Prototype reuse (original unsafe path)
print("Benchmarking prototype reuse + Setup + Energy...")
t0 = time.time()
for i in range(N):
    mol = mols[i % len(mols)]
    ff_proto.Setup(mol)
    ff_proto.Energy()
t1 = time.time()
rate_reuse = N / (t1 - t0)
print(f"  Prototype reuse: {N} runs in {t1-t0:.3f}s = {rate_reuse:.1f} mol/s")

overhead = rate_reuse / rate_newinstance
print(f"\nOverhead of MakeNewInstance vs reuse: {overhead:.2f}x")
print(f"MakeNewInstance achieves {rate_newinstance/rate_reuse*100:.1f}% of prototype reuse speed")
