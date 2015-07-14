# GenomeClustering

## Setup

* ```git clone``` the project
* ```cmake .``` to generate Makefile
* ```make``` to make executable
* ```./out/GenomeClustering <target fasta> <output path> <id> [<args>]``` to run the executable

## Usage

```
usage: ./out/GenomeClustering <fasta in> <cluster out> <similarity> [<args>]

(--cache_size | -c) <size>
    Set total cache size.
(--lru | -r) <size>
    Set the size of the LRU cache.
    Behaviour undefined with --cache_size.
(--lfu | -f) <size>
    Set the size of the LFU cache.
```

Example use:
```
./out/GenomeClustering res/p3_clean_C-148-2-Caecum_S128.fa clusters.uc 0.97
```

