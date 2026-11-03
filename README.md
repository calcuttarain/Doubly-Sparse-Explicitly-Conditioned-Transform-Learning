# Learning Doubly Sparse Explicitly Conditioned Transforms
This repository contains the MATLAB implementation for the _Learning Doubly Sparse Explicitly Conditioned Transforms_ paper.

## Code and Usage
- The proposed method is located in the `./TLAlgorithms/DoublySparseConditionedTL.m` file and its auxiliary routines are located in `./TLAlgorithms/DoublyConditionedTLRoutines/` folder.

- To test convergence, run `./test_convergence`.
- To test denoising, run `./test_denoising`.

## Citation
Please cite this work as:
```bibtex
@article{Pistol2026,
  title = {Learning Doubly Sparse Explicitly Conditioned Transforms},
  author = {Pistol, Tudor},
  journal = {Procedia Computer Science},
  volume = {292},
  pages = {7470--7479},
  year = {2026},
  doi = {10.1016/j.procs.2026.10.318}
}
```

## Acknowledgments
We thank the authors of the baseline frameworks and implementations for making their code publicly available.
