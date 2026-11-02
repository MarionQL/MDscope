"""Motif contact-query reuse microbenchmark (not whole-trajectory throughput).

Run: python scripts/benchmark_contact_reuse.py
80,000 lipid atoms, eight overlapping-region-sized 40-atom motifs, ten frames.
Reports Python allocation peak, not native process RSS. Checks exact result
checksums; timing is evidence, never a pass/fail threshold.
"""
import time,json,tracemalloc
import numpy as np
from MDAnalysis.lib.distances import capped_distance
rng=np.random.default_rng(2026)
box=np.array([160.,160.,160.,90.,90.,90.])
lipids=rng.uniform([0,0,35],[160,160,65],(80000,3)).astype(np.float32)
motifs=[rng.uniform([65,65,42],[95,95,58],(40,3)).astype(np.float32) for _ in range(8)]
results={};reference=None
for mode in ('before','after'):
 tracemalloc.start(); start=time.perf_counter();checksum=0
 for frame in range(10):
  cache=[]
  for xyz in motifs:
   if mode=='before':
    pairs=capped_distance(xyz,lipids,5.,box=box,return_distances=False)
   else:
    pairs,dist=capped_distance(xyz,lipids,5.,box=box,return_distances=True)
    cache.append((pairs,dist))
   checksum+=len(pairs)
  for i,xyz in enumerate(motifs):
   pairs,dist=(capped_distance(xyz,lipids,5.,box=box,return_distances=True)
               if mode=='before' else cache[i])
   checksum+=len(pairs)+float(dist.sum())
 elapsed=time.perf_counter()-start;_,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
 results[mode]={'seconds':elapsed,'python_peak_MiB':peak/2**20,'query_count':160 if mode=='before' else 80,'checksum':checksum}
 if reference is None:reference=checksum
 else:np.testing.assert_allclose(checksum,reference,rtol=0,atol=0)
results['speedup']=results['before']['seconds']/results['after']['seconds']
print(json.dumps(results,indent=2))
