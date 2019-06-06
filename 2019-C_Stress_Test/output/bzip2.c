#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <utime.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/times.h>
/*-------------------------------------------------------------*/
/*--- Public header file for the library.                   ---*/
/*---                                               bzlib.h ---*/
/*-------------------------------------------------------------*/
/*--
This file is a part of bzip2 and/or libbzip2, a program and
library for lossless, block-sorting data compression.

Copyright (C) 1996-2002 Julian R Seward.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software.  If you use this
software in a product, an acknowledgment in the product
documentation would be appreciated but is not required.

3. Altered source versions must be plainly marked as such, and must
not be misrepresented as being the original software.

4. The name of the author may not be used to endorse or promote
products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Julian Seward, Cambridge, UK.
jseward@acm.org
bzip2/libbzip2 version 1.0 of 21 March 2000

This program is based on (at least) the work of:
Mike Burrows
David Wheeler
Peter Fenwick
Alistair Moffat
Radford Neal
Ian H. Witten
Robert Sedgewick
Jon L. Bentley

For more information on these sources, see the manual.
--*/

struct anon_bzip2_c_88 {
   char *next_in;
   unsigned int avail_in;
   unsigned int total_in_lo32;
   unsigned int total_in_hi32;
   char *next_out;
   unsigned int avail_out;
   unsigned int total_out_lo32;
   unsigned int total_out_hi32;
   void *state;
   void *(*bzalloc) (void *, int, int);
   void (*bzfree) (void *, void *);
   void *opaque;
};

typedef struct anon_bzip2_c_88 bz_stream;
/*Need a definitition for FILE*/
/*windows.h define small to char*/
/*import windows dll dynamically*/
/*-- Core (low-level) library functions --*/
int clava_dcg_global[ 326 ] = {0};
extern int BZ2_bzCompressInit(bz_stream *strm, int blockSize100k, int verbosity, int workFactor);
extern int BZ2_bzCompress(bz_stream *strm, int action);
extern int BZ2_bzCompressEnd(bz_stream *strm);
extern int BZ2_bzDecompressInit(bz_stream *strm, int verbosity, int small);
extern int BZ2_bzDecompress(bz_stream *strm);
extern int BZ2_bzDecompressEnd(bz_stream *strm);
/*-- High(er) level library functions --*/

typedef void BZFILE;
extern BZFILE * BZ2_bzReadOpen(int *bzerror, FILE *f, int verbosity, int small, void *unused, int nUnused);
extern void BZ2_bzReadClose(int *bzerror, BZFILE *b);
extern void BZ2_bzReadGetUnused(int *bzerror, BZFILE *b, void **unused, int *nUnused);
extern int BZ2_bzRead(int *bzerror, BZFILE *b, void *buf, int len);
extern BZFILE * BZ2_bzWriteOpen(int *bzerror, FILE *f, int blockSize100k, int verbosity, int workFactor);
extern void BZ2_bzWrite(int *bzerror, BZFILE *b, void *buf, int len);
extern void BZ2_bzWriteClose(int *bzerror, BZFILE *b, int abandon, unsigned int *nbytes_in, unsigned int *nbytes_out);
extern void BZ2_bzWriteClose64(int *bzerror, BZFILE *b, int abandon, unsigned int *nbytes_in_lo32, unsigned int *nbytes_in_hi32, unsigned int *nbytes_out_lo32, unsigned int *nbytes_out_hi32);
/*-- Utility functions --*/
extern int BZ2_bzBuffToBuffCompress(char *dest, unsigned int *destLen, char *source, unsigned int sourceLen, int blockSize100k, int verbosity, int workFactor);
extern int BZ2_bzBuffToBuffDecompress(char *dest, unsigned int *destLen, char *source, unsigned int sourceLen, int small, int verbosity);
/*--
Code contributed by Yoshioka Tsuneo
(QWF00133@niftyserve.or.jp/tsuneo-y@is.aist-nara.ac.jp),
to support better zlib compatibility.
This code is not _officially_ part of libbzip2 (yet);
I haven't tested it, documented it, or considered the
threading-safeness of it.
If this code breaks, please contact both Yoshioka and me.
--*/
extern char const * BZ2_bzlibVersion();
extern BZFILE * BZ2_bzopen(char const *path, char const *mode);
extern BZFILE * BZ2_bzdopen(int fd, char const *mode);
extern int BZ2_bzread(BZFILE *b, void *buf, int len);
extern int BZ2_bzwrite(BZFILE *b, void *buf, int len);
extern int BZ2_bzflush(BZFILE *b);
extern void BZ2_bzclose(BZFILE *b);
extern char const * BZ2_bzerror(BZFILE *b, int *errnum);
/*-------------------------------------------------------------*/

/*--- end                                           bzlib.h ---*/

/*-------------------------------------------------------------*/

/*-------------------------------------------------------------*/

/*--- Private header file for the library.                  ---*/

/*---                                       bzlib_private.h ---*/

/*-------------------------------------------------------------*/

/*-- General stuff. --*/

typedef char Char;
typedef unsigned char Bool;
typedef unsigned char UChar;
typedef int Int32;
typedef unsigned int UInt32;
typedef short Int16;
typedef unsigned short UInt16;
/**/
extern void BZ2_bz__AssertH__fail(int errcode);
/**/

/**/

/**/

/**/

/**/

/**/

/**/

/**/

/*-- Header bytes. --*/

/*'B'*/

/*'Z'*/

/*'h'*/

/*'0'*/

/*-- Constants for the back end. --*/

/*-- Stuff for randomising repetitive blocks. --*/

extern Int32 BZ2_rNums[512];
/*-- Stuff for doing CRCs. --*/

extern UInt32 BZ2_crc32Table[256];
/*-- States and modes for compression. --*/
/*-- Structure holding all the compression-side stuff. --*/

struct anon_bzip2_c_492 {
   /*pointer back to the struct bz_stream*/
   bz_stream *strm;
   /*mode this stream is in, and whether inputting*/
   /*or outputting data*/
   Int32 mode;
   Int32 state;
   /*remembers avail_in when flush/finish requested*/
   UInt32 avail_in_expect;
   /*for doing the block sorting*/
   UInt32 *arr1;
   UInt32 *arr2;
   UInt32 *ftab;
   Int32 origPtr;
   /*aliases for arr1 and arr2*/
   UInt32 *ptr;
   UChar *block;
   UInt16 *mtfv;
   UChar *zbits;
   /*for deciding when to use the fallback sorting algorithm*/
   Int32 workFactor;
   /*run-length-encoding of the input*/
   UInt32 state_in_ch;
   Int32 state_in_len;
   Int32 rNToGo;
   Int32 rTPos;
   /*input and output limits and current posns*/
   Int32 nblock;
   Int32 nblockMAX;
   Int32 numZ;
   Int32 state_out_pos;
   /*map of bytes used in block*/
   Int32 nInUse;
   Bool inUse[256];
   UChar unseqToSeq[256];
   /*the buffer for bit stream creation*/
   UInt32 bsBuff;
   Int32 bsLive;
   /*block and combined CRCs*/
   UInt32 blockCRC;
   UInt32 combinedCRC;
   /*misc administratium*/
   Int32 verbosity;
   Int32 blockNo;
   Int32 blockSize100k;
   /*stuff for coding the MTF values*/
   Int32 nMTF;
   Int32 mtfFreq[258];
   UChar selector[18002];
   UChar selectorMtf[18002];
   UChar len[6][258];
   Int32 code[6][258];
   Int32 rfreq[6][258];
   /*second dimension: only 3 needed; 4 makes index calculations faster*/
   UInt32 len_pack[258][4];
};

typedef struct anon_bzip2_c_492 EState;
/*-- externs for compression. --*/
extern void BZ2_blockSort(EState *);
extern void BZ2_compressBlock(EState *, Bool);
extern void BZ2_bsInitWrite(EState *);
extern void BZ2_hbAssignCodes(Int32 *, UChar *, Int32, Int32, Int32);
extern void BZ2_hbMakeCodeLengths(UChar *, Int32 *, Int32, Int32);
/*-- states for decompression. --*/
/*-- Constants for the fast MTF decoder. --*/
/*-- Structure holding all the decompression-side stuff. --*/

struct anon_bzip2_c_643 {
   /*pointer back to the struct bz_stream*/
   bz_stream *strm;
   /*state indicator for this stream*/
   Int32 state;
   /*for doing the final run-length decoding*/
   UChar state_out_ch;
   Int32 state_out_len;
   Bool blockRandomised;
   Int32 rNToGo;
   Int32 rTPos;
   /*the buffer for bit stream reading*/
   UInt32 bsBuff;
   Int32 bsLive;
   /*misc administratium*/
   Int32 blockSize100k;
   Bool smallDecompress;
   Int32 currBlockNo;
   Int32 verbosity;
   /*for undoing the Burrows-Wheeler transform*/
   Int32 origPtr;
   UInt32 tPos;
   Int32 k0;
   Int32 unzftab[256];
   Int32 nblock_used;
   Int32 cftab[257];
   Int32 cftabCopy[257];
   /*for undoing the Burrows-Wheeler transform (FAST)*/
   UInt32 *tt;
   /*for undoing the Burrows-Wheeler transform (SMALL)*/
   UInt16 *ll16;
   UChar *ll4;
   /*stored and calculated CRCs*/
   UInt32 storedBlockCRC;
   UInt32 storedCombinedCRC;
   UInt32 calculatedBlockCRC;
   UInt32 calculatedCombinedCRC;
   /*map of bytes used in block*/
   Int32 nInUse;
   Bool inUse[256];
   Bool inUse16[16];
   UChar seqToUnseq[256];
   /*for decoding the MTF values*/
   UChar mtfa[4096];
   Int32 mtfbase[16];
   UChar selector[18002];
   UChar selectorMtf[18002];
   UChar len[6][258];
   Int32 limit[6][258];
   Int32 base[6][258];
   Int32 perm[6][258];
   Int32 minLens[6];
   /*save area for scalars in the main decompress code*/
   Int32 save_i;
   Int32 save_j;
   Int32 save_t;
   Int32 save_alphaSize;
   Int32 save_nGroups;
   Int32 save_nSelectors;
   Int32 save_EOB;
   Int32 save_groupNo;
   Int32 save_groupPos;
   Int32 save_nextSym;
   Int32 save_nblockMAX;
   Int32 save_nblock;
   Int32 save_es;
   Int32 save_N;
   Int32 save_curr;
   Int32 save_zt;
   Int32 save_zn;
   Int32 save_zvec;
   Int32 save_zj;
   Int32 save_gSel;
   Int32 save_gMinlen;
   Int32 *save_gLimit;
   Int32 *save_gBase;
   Int32 *save_gPerm;
};

typedef struct anon_bzip2_c_643 DState;
/*-- Macros for decompression. --*/
/*-- externs for decompression. --*/
extern Int32 BZ2_indexIntoF(Int32, Int32 *);
extern Int32 BZ2_decompress(DState *);
extern void BZ2_hbCreateDecodeTables(Int32 *, Int32 *, Int32 *, UChar *, Int32, Int32, Int32);
/*-- BZ_NO_STDIO seems to make NULL disappear on some platforms. --*/
/*-------------------------------------------------------------*/
/*--- end                                   bzlib_private.h ---*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*--- Block sorting machinery                               ---*/
/*---                                           blocksort.c ---*/
/*-------------------------------------------------------------*/
/*---------------------------------------------*/
/*--- Fallback O(N log(N)^2) sorting        ---*/
/*--- algorithm, for repetitive blocks      ---*/
/*---------------------------------------------*/
/*---------------------------------------------*/
inline static void fallbackSimpleSort(UInt32 *fmap, UInt32 *eclass, Int32 lo, Int32 hi) {
   Int32 i, j, tmp;
   UInt32 ec_tmp;
   if(lo == hi) 
   return;
   if(hi - lo > 3) {
      for(i = hi - 4; i >= lo; i--) {
         tmp = fmap[i];
         ec_tmp = eclass[tmp];
         for(j = i + 4; j <= hi && ec_tmp > eclass[fmap[j]]; j += 4) fmap[j - 4] = fmap[j];
         fmap[j - 4] = tmp;
      }
   }
   for(i = hi - 1; i >= lo; i--) {
      tmp = fmap[i];
      ec_tmp = eclass[tmp];
      for(j = i + 1; j <= hi && ec_tmp > eclass[fmap[j]]; j++) fmap[j - 1] = fmap[j];
      fmap[j - 1] = tmp;
   }
}

/*---------------------------------------------*/
static void fallbackQSort3(UInt32 *fmap, UInt32 *eclass, Int32 loSt, Int32 hiSt) {
   Int32 unLo, unHi, ltLo, gtHi, n, m;
   Int32 sp, lo, hi;
   UInt32 med, r, r3;
   Int32 stackLo[100];
   Int32 stackHi[100];
   r = 0;
   sp = 0;
   {
      stackLo[sp] = loSt;
      stackHi[sp] = hiSt;
      sp++;
   }
   ;
   while(sp > 0) {
      {
         if(!(sp < 100)) {
            clava_dcg_global[ 0 ]++;
            BZ2_bz__AssertH__fail(1004);
         }
      }
      ;
      {
         sp--;
         lo = stackLo[sp];
         hi = stackHi[sp];
      }
      ;
      if(hi - lo < 10) {
         clava_dcg_global[ 1 ]++;
         fallbackSimpleSort(fmap, eclass, lo, hi);
         continue;
      }
      /*Random partitioning.  Median of 3 sometimes fails to
      avoid bad cases.  Median of 9 seems to help but
      looks rather expensive.  This too seems to work but
      is cheaper.  Guidance for the magic constants
      7621 and 32768 is taken from Sedgewick's algorithms
      book, chapter 35.
      */
      r = ((r * 7621) + 1) % 32768;
      r3 = r % 3;
      if(r3 == 0) med = eclass[fmap[lo]];
      else if(r3 == 1) med = eclass[fmap[(lo + hi) >> 1]];
      else med = eclass[fmap[hi]];
      unLo = ltLo = lo;
      unHi = gtHi = hi;
      while(1) {
         while(1) {
            if(unLo > unHi) break;
            n = (Int32) eclass[fmap[unLo]] - (Int32) med;
            if(n == 0) {
               {
                  Int32 zztmp = fmap[unLo];
                  fmap[unLo] = fmap[ltLo];
                  fmap[ltLo] = zztmp;
               }
               ;
               ltLo++;
               unLo++;
               continue;
            }
            ;
            if(n > 0) break;
            unLo++;
         }
         while(1) {
            if(unLo > unHi) break;
            n = (Int32) eclass[fmap[unHi]] - (Int32) med;
            if(n == 0) {
               {
                  Int32 zztmp = fmap[unHi];
                  fmap[unHi] = fmap[gtHi];
                  fmap[gtHi] = zztmp;
               }
               ;
               gtHi--;
               unHi--;
               continue;
            }
            ;
            if(n < 0) break;
            unHi--;
         }
         if(unLo > unHi) break;
         {
            Int32 zztmp = fmap[unLo];
            fmap[unLo] = fmap[unHi];
            fmap[unHi] = zztmp;
         }
         ;
         unLo++;
         unHi--;
      }
      ;
      if(gtHi < ltLo) continue;
      n = ((ltLo - lo) < (unLo - ltLo)) ? (ltLo - lo) : (unLo - ltLo);
      {
         Int32 yyp1 = (lo);
         Int32 yyp2 = (unLo - n);
         Int32 yyn = (n);
         while(yyn > 0) {
            {
               Int32 zztmp = fmap[yyp1];
               fmap[yyp1] = fmap[yyp2];
               fmap[yyp2] = zztmp;
            }
            ;
            yyp1++;
            yyp2++;
            yyn--;
         }
      }
      ;
      m = ((hi - gtHi) < (gtHi - unHi)) ? (hi - gtHi) : (gtHi - unHi);
      {
         Int32 yyp1 = (unLo);
         Int32 yyp2 = (hi - m + 1);
         Int32 yyn = (m);
         while(yyn > 0) {
            {
               Int32 zztmp = fmap[yyp1];
               fmap[yyp1] = fmap[yyp2];
               fmap[yyp2] = zztmp;
            }
            ;
            yyp1++;
            yyp2++;
            yyn--;
         }
      }
      ;
      n = lo + unLo - ltLo - 1;
      m = hi - (gtHi - unHi) + 1;
      if(n - lo > hi - m) {
         {
            stackLo[sp] = lo;
            stackHi[sp] = n;
            sp++;
         }
         ;
         {
            stackLo[sp] = m;
            stackHi[sp] = hi;
            sp++;
         }
         ;
      }
      else {
         {
            stackLo[sp] = m;
            stackHi[sp] = hi;
            sp++;
         }
         ;
         {
            stackLo[sp] = lo;
            stackHi[sp] = n;
            sp++;
         }
         ;
      }
   }
}

/*---------------------------------------------*/
/*Pre:
nblock > 0
eclass exists for [0 .. nblock-1]
((UChar*)eclass) [0 .. nblock-1] holds block
ptr exists for [0 .. nblock-1]

Post:
((UChar*)eclass) [0 .. nblock-1] holds block
All other areas of eclass destroyed
fmap [0 .. nblock-1] holds sorted order
bhtab [ 0 .. 2+(nblock/32) ] destroyed
*/
static void fallbackSort(UInt32 *fmap, UInt32 *eclass, UInt32 *bhtab, Int32 nblock, Int32 verb) {
   Int32 ftab[257];
   Int32 ftabCopy[256];
   Int32 H, i, j, k, l, r, cc, cc1;
   Int32 nNotDone;
   Int32 nBhtab;
   UChar *eclass8 = (UChar *) eclass;
   /*--
   Initial 1-char radix sort to generate
   initial fmap and initial BH bits.
   --*/
   if(verb >= 4) {
      clava_dcg_global[ 2 ]++;
      fprintf(stderr, "        bucket sorting ...\n");
   }
   for(i = 0; i < 257; i++) ftab[i] = 0;
   for(i = 0; i < nblock; i++) ftab[eclass8[i]]++;
   for(i = 0; i < 256; i++) ftabCopy[i] = ftab[i];
   for(i = 1; i < 257; i++) ftab[i] += ftab[i - 1];
   for(i = 0; i < nblock; i++) {
      j = eclass8[i];
      k = ftab[j] - 1;
      ftab[j] = k;
      fmap[k] = i;
   }
   nBhtab = 2 + (nblock / 32);
   for(i = 0; i < nBhtab; i++) bhtab[i] = 0;
   for(i = 0; i < 256; i++) bhtab[(ftab[i]) >> 5] |= (1 << ((ftab[i]) & 31));
   /*--
   Inductively refine the buckets.  Kind-of an
   "exponential radix sort" (!), inspired by the
   Manber-Myers suffix array construction algorithm.
   --*/
   /*-- set sentinel bits for block-end detection --*/
   for(i = 0; i < 32; i++) {
      bhtab[(nblock + 2 * i) >> 5] |= (1 << ((nblock + 2 * i) & 31));
      bhtab[(nblock + 2 * i + 1) >> 5] &= ~(1 << ((nblock + 2 * i + 1) & 31));
   }
   /*-- the log(N) loop --*/
   H = 1;
   while(1) {
      if(verb >= 4) {
         clava_dcg_global[ 2 ]++;
         fprintf(stderr, "        depth %6d has ", H);
      }
      j = 0;
      for(i = 0; i < nblock; i++) {
         if((bhtab[(i) >> 5] & (1 << ((i) & 31)))) j = i;
         k = fmap[i] - H;
         if(k < 0) k += nblock;
         eclass[k] = j;
      }
      nNotDone = 0;
      r = -1;
      while(1) {
         /*-- find the next non-singleton bucket --*/
         k = r + 1;
         while((bhtab[(k) >> 5] & (1 << ((k) & 31))) && ((k) & 0x01f)) k++;
         if((bhtab[(k) >> 5] & (1 << ((k) & 31)))) {
            while(bhtab[(k) >> 5] == 0xffffffff) k += 32;
            while((bhtab[(k) >> 5] & (1 << ((k) & 31)))) k++;
         }
         l = k - 1;
         if(l >= nblock) break;
         while(!(bhtab[(k) >> 5] & (1 << ((k) & 31))) && ((k) & 0x01f)) k++;
         if(!(bhtab[(k) >> 5] & (1 << ((k) & 31)))) {
            while(bhtab[(k) >> 5] == 0x00000000) k += 32;
            while(!(bhtab[(k) >> 5] & (1 << ((k) & 31)))) k++;
         }
         r = k - 1;
         if(r >= nblock) break;
         /*-- now [l, r] bracket current bucket --*/
         if(r > l) {
            nNotDone += (r - l + 1);
            clava_dcg_global[ 3 ]++;
            fallbackQSort3(fmap, eclass, l, r);
            /*-- scan bucket and generate header bits--*/
            cc = -1;
            for(i = l; i <= r; i++) {
               cc1 = eclass[fmap[i]];
               if(cc != cc1) {
                  bhtab[(i) >> 5] |= (1 << ((i) & 31));
                  cc = cc1;
               }
               ;
            }
         }
      }
      if(verb >= 4) {
         clava_dcg_global[ 2 ]++;
         fprintf(stderr, "%6d unresolved strings\n", nNotDone);
      }
      H *= 2;
      if(H > nblock || nNotDone == 0) break;
   }
   /*--
   Reconstruct the original block in
   eclass8 [0 .. nblock-1], since the
   previous phase destroyed it.
   --*/
   if(verb >= 4) {
      clava_dcg_global[ 2 ]++;
      fprintf(stderr, "        reconstructing block ...\n");
   }
   j = 0;
   for(i = 0; i < nblock; i++) {
      while(ftabCopy[j] == 0) j++;
      ftabCopy[j]--;
      eclass8[fmap[i]] = (UChar) j;
   }
   {
      if(!(j < 256)) {
         clava_dcg_global[ 4 ]++;
         BZ2_bz__AssertH__fail(1005);
      }
   }
   ;
}

/*---------------------------------------------*/
/*--- The main, O(N^2 log(N)) sorting       ---*/
/*--- algorithm.  Faster for "normal"       ---*/
/*--- non-repetitive blocks.                ---*/
/*---------------------------------------------*/
/*---------------------------------------------*/
inline static Bool mainGtU(UInt32 i1, UInt32 i2, UChar *block, UInt16 *quadrant, UInt32 nblock, Int32 *budget) {
   Int32 k;
   UChar c1, c2;
   UInt16 s1, s2;
   ;
   /*1*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*2*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*3*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*4*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*5*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*6*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*7*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*8*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*9*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*10*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*11*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   /*12*/
   c1 = block[i1];
   c2 = block[i2];
   if(c1 != c2) 
   return (c1 > c2);
   i1++;
   i2++;
   k = nblock + 8;
   do  {
      /*1*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*2*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*3*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*4*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*5*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*6*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*7*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      /*8*/
      c1 = block[i1];
      c2 = block[i2];
      if(c1 != c2) 
      return (c1 > c2);
      s1 = quadrant[i1];
      s2 = quadrant[i2];
      if(s1 != s2) 
      return (s1 > s2);
      i1++;
      i2++;
      if(i1 >= nblock) i1 -= nblock;
      if(i2 >= nblock) i2 -= nblock;
      k -= 8;
      (*budget)--;
   }
   while (k >= 0);
   
   return ((Bool) 0);
}

/*---------------------------------------------*/

/*--
Knuth's increments seem to work better
than Incerpi-Sedgewick here.  Possibly
because the number of elems to sort is
usually small, typically <= 20.
--*/

static Int32 incs[14] = {1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 29524, 88573, 265720, 797161, 2391484};
static void mainSimpleSort(UInt32 *ptr, UChar *block, UInt16 *quadrant, Int32 nblock, Int32 lo, Int32 hi, Int32 d, Int32 *budget) {
   Int32 i, j, h, bigN, hp;
   UInt32 v;
   bigN = hi - lo + 1;
   if(bigN < 2) 
   return;
   hp = 0;
   while(incs[hp] < bigN) hp++;
   hp--;
   for(; hp >= 0; hp--) {
      h = incs[hp];
      i = lo + h;
      while(((Bool) 1)) {
         /*-- copy 1 --*/
         if(i > hi) break;
         v = ptr[i];
         j = i;
         while(mainGtU(ptr[j - h] + d, v + d, block, quadrant, nblock, budget)) {
            clava_dcg_global[ 5 ]++;
            ptr[j] = ptr[j - h];
            j = j - h;
            if(j <= (lo + h - 1)) break;
         }
         ptr[j] = v;
         i++;
         /*-- copy 2 --*/
         if(i > hi) break;
         v = ptr[i];
         j = i;
         while(mainGtU(ptr[j - h] + d, v + d, block, quadrant, nblock, budget)) {
            clava_dcg_global[ 5 ]++;
            ptr[j] = ptr[j - h];
            j = j - h;
            if(j <= (lo + h - 1)) break;
         }
         ptr[j] = v;
         i++;
         /*-- copy 3 --*/
         if(i > hi) break;
         v = ptr[i];
         j = i;
         while(mainGtU(ptr[j - h] + d, v + d, block, quadrant, nblock, budget)) {
            clava_dcg_global[ 5 ]++;
            ptr[j] = ptr[j - h];
            j = j - h;
            if(j <= (lo + h - 1)) break;
         }
         ptr[j] = v;
         i++;
         if(*budget < 0) 
         return;
      }
   }
}

/*---------------------------------------------*/
/*--
The following is an implementation of
an elegant 3-way quicksort for strings,
described in a paper "Fast Algorithms for
Sorting and Searching Strings", by Robert
Sedgewick and Jon L. Bentley.
--*/
inline static UChar mmed3(UChar a, UChar b, UChar c) {
   UChar t;
   if(a > b) {
      t = a;
      a = b;
      b = t;
   }
   ;
   if(b > c) {
      b = c;
      if(a > b) b = a;
   }
   
   return b;
}

static void mainQSort3(UInt32 *ptr, UChar *block, UInt16 *quadrant, Int32 nblock, Int32 loSt, Int32 hiSt, Int32 dSt, Int32 *budget) {
   Int32 unLo, unHi, ltLo, gtHi, n, m, med;
   Int32 sp, lo, hi, d;
   Int32 stackLo[100];
   Int32 stackHi[100];
   Int32 stackD[100];
   Int32 nextLo[3];
   Int32 nextHi[3];
   Int32 nextD[3];
   sp = 0;
   {
      stackLo[sp] = loSt;
      stackHi[sp] = hiSt;
      stackD[sp] = dSt;
      sp++;
   }
   ;
   while(sp > 0) {
      {
         if(!(sp < 100)) {
            clava_dcg_global[ 6 ]++;
            BZ2_bz__AssertH__fail(1001);
         }
      }
      ;
      {
         sp--;
         lo = stackLo[sp];
         hi = stackHi[sp];
         d = stackD[sp];
      }
      ;
      if(hi - lo < 20 || d > (2 + 12)) {
         clava_dcg_global[ 7 ]++;
         mainSimpleSort(ptr, block, quadrant, nblock, lo, hi, d, budget);
         if(*budget < 0) 
         return;
         continue;
      }
      clava_dcg_global[ 8 ]++;
      med = (Int32) mmed3(block[ptr[lo] + d], block[ptr[hi] + d], block[ptr[(lo + hi) >> 1] + d]);
      unLo = ltLo = lo;
      unHi = gtHi = hi;
      while(((Bool) 1)) {
         while(((Bool) 1)) {
            if(unLo > unHi) break;
            n = ((Int32) block[ptr[unLo] + d]) - med;
            if(n == 0) {
               {
                  Int32 zztmp = ptr[unLo];
                  ptr[unLo] = ptr[ltLo];
                  ptr[ltLo] = zztmp;
               }
               ;
               ltLo++;
               unLo++;
               continue;
            }
            ;
            if(n > 0) break;
            unLo++;
         }
         while(((Bool) 1)) {
            if(unLo > unHi) break;
            n = ((Int32) block[ptr[unHi] + d]) - med;
            if(n == 0) {
               {
                  Int32 zztmp = ptr[unHi];
                  ptr[unHi] = ptr[gtHi];
                  ptr[gtHi] = zztmp;
               }
               ;
               gtHi--;
               unHi--;
               continue;
            }
            ;
            if(n < 0) break;
            unHi--;
         }
         if(unLo > unHi) break;
         {
            Int32 zztmp = ptr[unLo];
            ptr[unLo] = ptr[unHi];
            ptr[unHi] = zztmp;
         }
         ;
         unLo++;
         unHi--;
      }
      ;
      if(gtHi < ltLo) {
         {
            stackLo[sp] = lo;
            stackHi[sp] = hi;
            stackD[sp] = d + 1;
            sp++;
         }
         ;
         continue;
      }
      n = ((ltLo - lo) < (unLo - ltLo)) ? (ltLo - lo) : (unLo - ltLo);
      {
         Int32 yyp1 = (lo);
         Int32 yyp2 = (unLo - n);
         Int32 yyn = (n);
         while(yyn > 0) {
            {
               Int32 zztmp = ptr[yyp1];
               ptr[yyp1] = ptr[yyp2];
               ptr[yyp2] = zztmp;
            }
            ;
            yyp1++;
            yyp2++;
            yyn--;
         }
      }
      ;
      m = ((hi - gtHi) < (gtHi - unHi)) ? (hi - gtHi) : (gtHi - unHi);
      {
         Int32 yyp1 = (unLo);
         Int32 yyp2 = (hi - m + 1);
         Int32 yyn = (m);
         while(yyn > 0) {
            {
               Int32 zztmp = ptr[yyp1];
               ptr[yyp1] = ptr[yyp2];
               ptr[yyp2] = zztmp;
            }
            ;
            yyp1++;
            yyp2++;
            yyn--;
         }
      }
      ;
      n = lo + unLo - ltLo - 1;
      m = hi - (gtHi - unHi) + 1;
      nextLo[0] = lo;
      nextHi[0] = n;
      nextD[0] = d;
      nextLo[1] = m;
      nextHi[1] = hi;
      nextD[1] = d;
      nextLo[2] = n + 1;
      nextHi[2] = m - 1;
      nextD[2] = d + 1;
      if((nextHi[0] - nextLo[0]) < (nextHi[1] - nextLo[1])) {
         Int32 tz;
         tz = nextLo[0];
         nextLo[0] = nextLo[1];
         nextLo[1] = tz;
         tz = nextHi[0];
         nextHi[0] = nextHi[1];
         nextHi[1] = tz;
         tz = nextD[0];
         nextD[0] = nextD[1];
         nextD[1] = tz;
      }
      ;
      if((nextHi[1] - nextLo[1]) < (nextHi[2] - nextLo[2])) {
         Int32 tz;
         tz = nextLo[1];
         nextLo[1] = nextLo[2];
         nextLo[2] = tz;
         tz = nextHi[1];
         nextHi[1] = nextHi[2];
         nextHi[2] = tz;
         tz = nextD[1];
         nextD[1] = nextD[2];
         nextD[2] = tz;
      }
      ;
      if((nextHi[0] - nextLo[0]) < (nextHi[1] - nextLo[1])) {
         Int32 tz;
         tz = nextLo[0];
         nextLo[0] = nextLo[1];
         nextLo[1] = tz;
         tz = nextHi[0];
         nextHi[0] = nextHi[1];
         nextHi[1] = tz;
         tz = nextD[0];
         nextD[0] = nextD[1];
         nextD[1] = tz;
      }
      ;
      ;
      ;
      {
         stackLo[sp] = nextLo[0];
         stackHi[sp] = nextHi[0];
         stackD[sp] = nextD[0];
         sp++;
      }
      ;
      {
         stackLo[sp] = nextLo[1];
         stackHi[sp] = nextHi[1];
         stackD[sp] = nextD[1];
         sp++;
      }
      ;
      {
         stackLo[sp] = nextLo[2];
         stackHi[sp] = nextHi[2];
         stackD[sp] = nextD[2];
         sp++;
      }
      ;
   }
}

/*---------------------------------------------*/
/*Pre:
nblock > N_OVERSHOOT
block32 exists for [0 .. nblock-1 +N_OVERSHOOT]
((UChar*)block32) [0 .. nblock-1] holds block
ptr exists for [0 .. nblock-1]

Post:
((UChar*)block32) [0 .. nblock-1] holds block
All other areas of block32 destroyed
ftab [0 .. 65536 ] destroyed
ptr [0 .. nblock-1] holds sorted order
if (*budget < 0), sorting was abandoned
*/
static void mainSort(UInt32 *ptr, UChar *block, UInt16 *quadrant, UInt32 *ftab, Int32 nblock, Int32 verb, Int32 *budget) {
   Int32 i, j, k, ss, sb;
   Int32 runningOrder[256];
   Bool bigDone[256];
   Int32 copyStart[256];
   Int32 copyEnd[256];
   UChar c1;
   Int32 numQSorted;
   UInt16 s;
   if(verb >= 4) {
      clava_dcg_global[ 9 ]++;
      fprintf(stderr, "        main sort initialise ...\n");
   }
   /*-- set up the 2-byte frequency table --*/
   for(i = 65536; i >= 0; i--) ftab[i] = 0;
   j = block[0] << 8;
   i = nblock - 1;
   for(; i >= 3; i -= 4) {
      quadrant[i] = 0;
      j = (j >> 8) | (((UInt16) block[i]) << 8);
      ftab[j]++;
      quadrant[i - 1] = 0;
      j = (j >> 8) | (((UInt16) block[i - 1]) << 8);
      ftab[j]++;
      quadrant[i - 2] = 0;
      j = (j >> 8) | (((UInt16) block[i - 2]) << 8);
      ftab[j]++;
      quadrant[i - 3] = 0;
      j = (j >> 8) | (((UInt16) block[i - 3]) << 8);
      ftab[j]++;
   }
   for(; i >= 0; i--) {
      quadrant[i] = 0;
      j = (j >> 8) | (((UInt16) block[i]) << 8);
      ftab[j]++;
   }
   /*-- (emphasises close relationship of block & quadrant) --*/
   for(i = 0; i < (2 + 12 + 18 + 2); i++) {
      block[nblock + i] = block[i];
      quadrant[nblock + i] = 0;
   }
   if(verb >= 4) {
      clava_dcg_global[ 9 ]++;
      fprintf(stderr, "        bucket sorting ...\n");
   }
   /*-- Complete the initial radix sort --*/
   for(i = 1; i <= 65536; i++) ftab[i] += ftab[i - 1];
   s = block[0] << 8;
   i = nblock - 1;
   for(; i >= 3; i -= 4) {
      s = (s >> 8) | (block[i] << 8);
      j = ftab[s] - 1;
      ftab[s] = j;
      ptr[j] = i;
      s = (s >> 8) | (block[i - 1] << 8);
      j = ftab[s] - 1;
      ftab[s] = j;
      ptr[j] = i - 1;
      s = (s >> 8) | (block[i - 2] << 8);
      j = ftab[s] - 1;
      ftab[s] = j;
      ptr[j] = i - 2;
      s = (s >> 8) | (block[i - 3] << 8);
      j = ftab[s] - 1;
      ftab[s] = j;
      ptr[j] = i - 3;
   }
   for(; i >= 0; i--) {
      s = (s >> 8) | (block[i] << 8);
      j = ftab[s] - 1;
      ftab[s] = j;
      ptr[j] = i;
   }
   /*--
   Now ftab contains the first loc of every small bucket.
   Calculate the running order, from smallest to largest
   big bucket.
   --*/
   for(i = 0; i <= 255; i++) {
      bigDone[i] = ((Bool) 0);
      runningOrder[i] = i;
   }
   {
      Int32 vv;
      Int32 h = 1;
      do  h = 3 * h + 1;while (h <= 256);
      do  {
         h = h / 3;
         for(i = h; i <= 255; i++) {
            vv = runningOrder[i];
            j = i;
            while((ftab[((runningOrder[j - h]) + 1) << 8] - ftab[(runningOrder[j - h]) << 8]) > (ftab[((vv) + 1) << 8] - ftab[(vv) << 8])) {
               runningOrder[j] = runningOrder[j - h];
               j = j - h;
               if(j <= (h - 1)) goto zero;
            }
            zero:
            runningOrder[j] = vv;
         }
      }
      while (h != 1);
   }
   /*--
   The main sorting loop.
   --*/
   numQSorted = 0;
   for(i = 0; i <= 255; i++) {
      /*--
      Process big buckets, starting with the least full.
      Basically this is a 3-step process in which we call
      mainQSort3 to sort the small buckets [ss, j], but
      also make a big effort to avoid the calls if we can.
      --*/
      ss = runningOrder[i];
      /*--
      Step 1:
      Complete the big bucket [ss] by quicksorting
      any unsorted small buckets [ss, j], for j != ss.
      Hopefully previous pointer-scanning phases have already
      completed many of the small buckets [ss, j], so
      we don't have to sort them at all.
      --*/
      for(j = 0; j <= 255; j++) {
         if(j != ss) {
            sb = (ss << 8) + j;
            if(!(ftab[sb] & (1 << 21))) {
               Int32 lo = ftab[sb] & (~((1 << 21)));
               Int32 hi = (ftab[sb + 1] & (~((1 << 21)))) - 1;
               if(hi > lo) {
                  if(verb >= 4) {
                     clava_dcg_global[ 9 ]++;
                     fprintf(stderr, "        qsort [0x%x, 0x%x]   done %d   this %d\n", ss, j, numQSorted, hi - lo + 1);
                  }
                  clava_dcg_global[ 10 ]++;
                  mainQSort3(ptr, block, quadrant, nblock, lo, hi, 2, budget);
                  numQSorted += (hi - lo + 1);
                  if(*budget < 0) 
                  return;
               }
            }
            ftab[sb] |= (1 << 21);
         }
      }
      {
         if(!(!bigDone[ss])) {
            clava_dcg_global[ 11 ]++;
            BZ2_bz__AssertH__fail(1006);
         }
      }
      ;
      /*--
      Step 2:
      Now scan this big bucket [ss] so as to synthesise the
      sorted order for small buckets [t, ss] for all t,
      including, magically, the bucket [ss,ss] too.
      This will avoid doing Real Work in subsequent Step 1's.
      --*/
      {
         for(j = 0; j <= 255; j++) {
            copyStart[j] = ftab[(j << 8) + ss] & (~((1 << 21)));
            copyEnd[j] = (ftab[(j << 8) + ss + 1] & (~((1 << 21)))) - 1;
         }
         for(j = ftab[ss << 8] & (~((1 << 21))); j < copyStart[ss]; j++) {
            k = ptr[j] - 1;
            if(k < 0) k += nblock;
            c1 = block[k];
            if(!bigDone[c1]) ptr[copyStart[c1]++] = k;
         }
         for(j = (ftab[(ss + 1) << 8] & (~((1 << 21)))) - 1; j > copyEnd[ss]; j--) {
            k = ptr[j] - 1;
            if(k < 0) k += nblock;
            c1 = block[k];
            if(!bigDone[c1]) ptr[copyEnd[c1]--] = k;
         }
      }
      {
         if(!((copyStart[ss] - 1 == copyEnd[ss]) || (copyStart[ss] == 0 && copyEnd[ss] == nblock - 1))) {
            clava_dcg_global[ 11 ]++;
            BZ2_bz__AssertH__fail(1007);
         }
      }
      /*Extremely rare case missing in bzip2-1.0.0 and 1.0.1.
      Necessity for this case is demonstrated by compressing
      a sequence of approximately 48.5 million of character
      251; 1.0.0/1.0.1 will then die here.*/
      for(j = 0; j <= 255; j++) ftab[(j << 8) + ss] |= (1 << 21);
      /*--
      Step 3:
      The [ss] big bucket is now done.  Record this fact,
      and update the quadrant descriptors.  Remember to
      update quadrants in the overshoot area too, if
      necessary.  The "if (i < 255)" test merely skips
      this updating for the last bucket processed, since
      updating for the last bucket is pointless.
      
      The quadrant array provides a way to incrementally
      cache sort orderings, as they appear, so as to
      make subsequent comparisons in fullGtU() complete
      faster.  For repetitive blocks this makes a big
      difference (but not big enough to be able to avoid
      the fallback sorting mechanism, exponential radix sort).
      
      The precise meaning is: at all times:
      
      for 0 <= i < nblock and 0 <= j <= nblock
      
      if block[i] != block[j],
      
      then the relative values of quadrant[i] and
      quadrant[j] are meaningless.
      
      else {
      if quadrant[i] < quadrant[j]
      then the string starting at i lexicographically
      precedes the string starting at j
      
      else if quadrant[i] > quadrant[j]
      then the string starting at j lexicographically
      precedes the string starting at i
      
      else
      the relative ordering of the strings starting
      at i and j has not yet been determined.
      }
      --*/
      bigDone[ss] = ((Bool) 1);
      if(i < 255) {
         Int32 bbStart = ftab[ss << 8] & (~((1 << 21)));
         Int32 bbSize = (ftab[(ss + 1) << 8] & (~((1 << 21)))) - bbStart;
         Int32 shifts = 0;
         while((bbSize >> shifts) > 65534) shifts++;
         for(j = bbSize - 1; j >= 0; j--) {
            Int32 a2update = ptr[bbStart + j];
            UInt16 qVal = (UInt16) (j >> shifts);
            quadrant[a2update] = qVal;
            if(a2update < (2 + 12 + 18 + 2)) quadrant[a2update + nblock] = qVal;
         }
         {
            if(!(((bbSize - 1) >> shifts) <= 65535)) {
               clava_dcg_global[ 11 ]++;
               BZ2_bz__AssertH__fail(1002);
            }
         }
         ;
      }
   }
   if(verb >= 4) {
      clava_dcg_global[ 9 ]++;
      fprintf(stderr, "        %d pointers, %d sorted, %d scanned\n", nblock, numQSorted, nblock - numQSorted);
   }
}

/*---------------------------------------------*/
/*Pre:
nblock > 0
arr2 exists for [0 .. nblock-1 +N_OVERSHOOT]
((UChar*)arr2)  [0 .. nblock-1] holds block
arr1 exists for [0 .. nblock-1]

Post:
((UChar*)arr2) [0 .. nblock-1] holds block
All other areas of block destroyed
ftab [ 0 .. 65536 ] destroyed
arr1 [0 .. nblock-1] holds sorted order
*/
void BZ2_blockSort(EState *s) {
   UInt32 *ptr = s->ptr;
   UChar *block = s->block;
   UInt32 *ftab = s->ftab;
   Int32 nblock = s->nblock;
   Int32 verb = s->verbosity;
   Int32 wfact = s->workFactor;
   UInt16 *quadrant;
   Int32 budget;
   Int32 budgetInit;
   Int32 i;
   if(nblock < 10000) {
      clava_dcg_global[ 12 ]++;
      fallbackSort(s->arr1, s->arr2, ftab, nblock, verb);
   }
   else {
      /*Calculate the location for quadrant, remembering to get
      the alignment right.  Assumes that &(block[0]) is at least
      2-byte aligned -- this should be ok since block is really
      the first section of arr2.
      */
      i = nblock + (2 + 12 + 18 + 2);
      if(i & 1) i++;
      quadrant = (UInt16 *) (&(block[i]));
      /*(wfact-1) / 3 puts the default-factor-30
      transition point at very roughly the same place as
      with v0.1 and v0.9.0.
      Not that it particularly matters any more, since the
      resulting compressed stream is now the same regardless
      of whether or not we use the main sort or fallback sort.
      */
      if(wfact < 1) wfact = 1;
      if(wfact > 100) wfact = 100;
      budgetInit = nblock * ((wfact - 1) / 3);
      budget = budgetInit;
      clava_dcg_global[ 13 ]++;
      mainSort(ptr, block, quadrant, ftab, nblock, verb, &budget);
      if(verb >= 3) {
         clava_dcg_global[ 14 ]++;
         fprintf(stderr, "      %d work, %d block, ratio %5.2f\n", budgetInit - budget, nblock, (float) (budgetInit - budget) / (float) (nblock == 0 ? 1 : nblock));
      }
      if(budget < 0) {
         if(verb >= 2) {
            clava_dcg_global[ 14 ]++;
            fprintf(stderr, "    too repetitive; using fallback sorting algorithm\n");
         }
         clava_dcg_global[ 12 ]++;
         fallbackSort(s->arr1, s->arr2, ftab, nblock, verb);
      }
   }
   s->origPtr = -1;
   for(i = 0; i < s->nblock; i++) if(ptr[i] == 0) {
      s->origPtr = i;
      break;
   }
   ;
   {
      if(!(s->origPtr != -1)) {
         clava_dcg_global[ 15 ]++;
         BZ2_bz__AssertH__fail(1003);
      }
   }
   ;
}

/*-------------------------------------------------------------*/
/*--- end                                       blocksort.c ---*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*--- Huffman coding low-level stuff                        ---*/
/*---                                             huffman.c ---*/
/*-------------------------------------------------------------*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
void BZ2_hbMakeCodeLengths(UChar *len, Int32 *freq, Int32 alphaSize, Int32 maxLen) {
   /*--
   Nodes and heap entries run from 1.  Entry 0
   for both the heap and nodes is a sentinel.
   --*/
   Int32 nNodes, nHeap, n1, n2, i, j, k;
   Bool tooLong;
   Int32 heap[260];
   Int32 weight[516];
   Int32 parent[516];
   for(i = 0; i < alphaSize; i++) weight[i + 1] = (freq[i] == 0 ? 1 : freq[i]) << 8;
   while(((Bool) 1)) {
      nNodes = alphaSize;
      nHeap = 0;
      heap[0] = 0;
      weight[0] = 0;
      parent[0] = -2;
      for(i = 1; i <= alphaSize; i++) {
         parent[i] = -1;
         nHeap++;
         heap[nHeap] = i;
         {
            Int32 zz, tmp;
            zz = nHeap;
            tmp = heap[zz];
            while(weight[tmp] < weight[heap[zz >> 1]]) {
               heap[zz] = heap[zz >> 1];
               zz >>= 1;
            }
            heap[zz] = tmp;
         }
         ;
      }
      {
         if(!(nHeap < (258 + 2))) {
            clava_dcg_global[ 16 ]++;
            BZ2_bz__AssertH__fail(2001);
         }
      }
      ;
      while(nHeap > 1) {
         n1 = heap[1];
         heap[1] = heap[nHeap];
         nHeap--;
         {
            Int32 zz, yy, tmp;
            zz = 1;
            tmp = heap[zz];
            while(((Bool) 1)) {
               yy = zz << 1;
               if(yy > nHeap) break;
               if(yy < nHeap && weight[heap[yy + 1]] < weight[heap[yy]]) yy++;
               if(weight[tmp] < weight[heap[yy]]) break;
               heap[zz] = heap[yy];
               zz = yy;
            }
            heap[zz] = tmp;
         }
         ;
         n2 = heap[1];
         heap[1] = heap[nHeap];
         nHeap--;
         {
            Int32 zz, yy, tmp;
            zz = 1;
            tmp = heap[zz];
            while(((Bool) 1)) {
               yy = zz << 1;
               if(yy > nHeap) break;
               if(yy < nHeap && weight[heap[yy + 1]] < weight[heap[yy]]) yy++;
               if(weight[tmp] < weight[heap[yy]]) break;
               heap[zz] = heap[yy];
               zz = yy;
            }
            heap[zz] = tmp;
         }
         ;
         nNodes++;
         parent[n1] = parent[n2] = nNodes;
         weight[nNodes] = (((weight[n1]) & 0xffffff00) + ((weight[n2]) & 0xffffff00)) | (1 + ((((weight[n1]) & 0x000000ff)) > (((weight[n2]) & 0x000000ff)) ? (((weight[n1]) & 0x000000ff)) : (((weight[n2]) & 0x000000ff))));
         parent[nNodes] = -1;
         nHeap++;
         heap[nHeap] = nNodes;
         {
            Int32 zz, tmp;
            zz = nHeap;
            tmp = heap[zz];
            while(weight[tmp] < weight[heap[zz >> 1]]) {
               heap[zz] = heap[zz >> 1];
               zz >>= 1;
            }
            heap[zz] = tmp;
         }
         ;
      }
      {
         if(!(nNodes < (258 * 2))) {
            clava_dcg_global[ 16 ]++;
            BZ2_bz__AssertH__fail(2002);
         }
      }
      ;
      tooLong = ((Bool) 0);
      for(i = 1; i <= alphaSize; i++) {
         j = 0;
         k = i;
         while(parent[k] >= 0) {
            k = parent[k];
            j++;
         }
         len[i - 1] = j;
         if(j > maxLen) tooLong = ((Bool) 1);
      }
      if(!tooLong) break;
      for(i = 1; i < alphaSize; i++) {
         j = weight[i] >> 8;
         j = 1 + (j / 2);
         weight[i] = j << 8;
      }
   }
}

/*---------------------------------------------------*/
void BZ2_hbAssignCodes(Int32 *code, UChar *length, Int32 minLen, Int32 maxLen, Int32 alphaSize) {
   Int32 n, vec, i;
   vec = 0;
   for(n = minLen; n <= maxLen; n++) {
      for(i = 0; i < alphaSize; i++) if(length[i] == n) {
         code[i] = vec;
         vec++;
      }
      ;
      vec <<= 1;
   }
}

/*---------------------------------------------------*/
void BZ2_hbCreateDecodeTables(Int32 *limit, Int32 *base, Int32 *perm, UChar *length, Int32 minLen, Int32 maxLen, Int32 alphaSize) {
   Int32 pp, i, j, vec;
   pp = 0;
   for(i = minLen; i <= maxLen; i++) for(j = 0; j < alphaSize; j++) if(length[j] == i) {
      perm[pp] = j;
      pp++;
   }
   ;
   for(i = 0; i < 23; i++) base[i] = 0;
   for(i = 0; i < alphaSize; i++) base[length[i] + 1]++;
   for(i = 1; i < 23; i++) base[i] += base[i - 1];
   for(i = 0; i < 23; i++) limit[i] = 0;
   vec = 0;
   for(i = minLen; i <= maxLen; i++) {
      vec += (base[i + 1] - base[i]);
      limit[i] = vec - 1;
      vec <<= 1;
   }
   for(i = minLen + 1; i <= maxLen; i++) base[i] = ((limit[i - 1] + 1) << 1) - base[i];
}

/*-------------------------------------------------------------*/

/*--- end                                         huffman.c ---*/

/*-------------------------------------------------------------*/

/*-------------------------------------------------------------*/

/*--- Table for doing CRCs                                  ---*/

/*---                                            crctable.c ---*/

/*-------------------------------------------------------------*/

/*--
I think this is an implementation of the AUTODIN-II,
Ethernet & FDDI 32-bit CRC standard.  Vaguely derived
from code by Rob Warnock, in Section 51 of the
comp.compression FAQ.
--*/

UInt32 BZ2_crc32Table[256] = {0x00000000L, 0x04c11db7L, 0x09823b6eL, 0x0d4326d9L, 0x130476dcL, 0x17c56b6bL, 0x1a864db2L, 0x1e475005L, 0x2608edb8L, 0x22c9f00fL, 0x2f8ad6d6L, 0x2b4bcb61L, 0x350c9b64L, 0x31cd86d3L, 0x3c8ea00aL, 0x384fbdbdL, 0x4c11db70L, 0x48d0c6c7L, 0x4593e01eL, 0x4152fda9L, 0x5f15adacL, 0x5bd4b01bL, 0x569796c2L, 0x52568b75L, 0x6a1936c8L, 0x6ed82b7fL, 0x639b0da6L, 0x675a1011L, 0x791d4014L, 0x7ddc5da3L, 0x709f7b7aL, 0x745e66cdL, 0x9823b6e0L, 0x9ce2ab57L, 0x91a18d8eL, 0x95609039L, 0x8b27c03cL, 0x8fe6dd8bL, 0x82a5fb52L, 0x8664e6e5L, 0xbe2b5b58L, 0xbaea46efL, 0xb7a96036L, 0xb3687d81L, 0xad2f2d84L, 0xa9ee3033L, 0xa4ad16eaL, 0xa06c0b5dL, 0xd4326d90L, 0xd0f37027L, 0xddb056feL, 0xd9714b49L, 0xc7361b4cL, 0xc3f706fbL, 0xceb42022L, 0xca753d95L, 0xf23a8028L, 0xf6fb9d9fL, 0xfbb8bb46L, 0xff79a6f1L, 0xe13ef6f4L, 0xe5ffeb43L, 0xe8bccd9aL, 0xec7dd02dL, 0x34867077L, 0x30476dc0L, 0x3d044b19L, 0x39c556aeL, 0x278206abL, 0x23431b1cL, 0x2e003dc5L, 0x2ac12072L, 0x128e9dcfL, 0x164f8078L, 0x1b0ca6a1L, 0x1fcdbb16L, 0x018aeb13L, 0x054bf6a4L, 0x0808d07dL, 0x0cc9cdcaL, 0x7897ab07L, 0x7c56b6b0L, 0x71159069L, 0x75d48ddeL, 0x6b93dddbL, 0x6f52c06cL, 0x6211e6b5L, 0x66d0fb02L, 0x5e9f46bfL, 0x5a5e5b08L, 0x571d7dd1L, 0x53dc6066L, 0x4d9b3063L, 0x495a2dd4L, 0x44190b0dL, 0x40d816baL, 0xaca5c697L, 0xa864db20L, 0xa527fdf9L, 0xa1e6e04eL, 0xbfa1b04bL, 0xbb60adfcL, 0xb6238b25L, 0xb2e29692L, 0x8aad2b2fL, 0x8e6c3698L, 0x832f1041L, 0x87ee0df6L, 0x99a95df3L, 0x9d684044L, 0x902b669dL, 0x94ea7b2aL, 0xe0b41de7L, 0xe4750050L, 0xe9362689L, 0xedf73b3eL, 0xf3b06b3bL, 0xf771768cL, 0xfa325055L, 0xfef34de2L, 0xc6bcf05fL, 0xc27dede8L, 0xcf3ecb31L, 0xcbffd686L, 0xd5b88683L, 0xd1799b34L, 0xdc3abdedL, 0xd8fba05aL, 0x690ce0eeL, 0x6dcdfd59L, 0x608edb80L, 0x644fc637L, 0x7a089632L, 0x7ec98b85L, 0x738aad5cL, 0x774bb0ebL, 0x4f040d56L, 0x4bc510e1L, 0x46863638L, 0x42472b8fL, 0x5c007b8aL, 0x58c1663dL, 0x558240e4L, 0x51435d53L, 0x251d3b9eL, 0x21dc2629L, 0x2c9f00f0L, 0x285e1d47L, 0x36194d42L, 0x32d850f5L, 0x3f9b762cL, 0x3b5a6b9bL, 0x0315d626L, 0x07d4cb91L, 0x0a97ed48L, 0x0e56f0ffL, 0x1011a0faL, 0x14d0bd4dL, 0x19939b94L, 0x1d528623L, 0xf12f560eL, 0xf5ee4bb9L, 0xf8ad6d60L, 0xfc6c70d7L, 0xe22b20d2L, 0xe6ea3d65L, 0xeba91bbcL, 0xef68060bL, 0xd727bbb6L, 0xd3e6a601L, 0xdea580d8L, 0xda649d6fL, 0xc423cd6aL, 0xc0e2d0ddL, 0xcda1f604L, 0xc960ebb3L, 0xbd3e8d7eL, 0xb9ff90c9L, 0xb4bcb610L, 0xb07daba7L, 0xae3afba2L, 0xaafbe615L, 0xa7b8c0ccL, 0xa379dd7bL, 0x9b3660c6L, 0x9ff77d71L, 0x92b45ba8L, 0x9675461fL, 0x8832161aL, 0x8cf30badL, 0x81b02d74L, 0x857130c3L, 0x5d8a9099L, 0x594b8d2eL, 0x5408abf7L, 0x50c9b640L, 0x4e8ee645L, 0x4a4ffbf2L, 0x470cdd2bL, 0x43cdc09cL, 0x7b827d21L, 0x7f436096L, 0x7200464fL, 0x76c15bf8L, 0x68860bfdL, 0x6c47164aL, 0x61043093L, 0x65c52d24L, 0x119b4be9L, 0x155a565eL, 0x18197087L, 0x1cd86d30L, 0x029f3d35L, 0x065e2082L, 0x0b1d065bL, 0x0fdc1becL, 0x3793a651L, 0x3352bbe6L, 0x3e119d3fL, 0x3ad08088L, 0x2497d08dL, 0x2056cd3aL, 0x2d15ebe3L, 0x29d4f654L, 0xc5a92679L, 0xc1683bceL, 0xcc2b1d17L, 0xc8ea00a0L, 0xd6ad50a5L, 0xd26c4d12L, 0xdf2f6bcbL, 0xdbee767cL, 0xe3a1cbc1L, 0xe760d676L, 0xea23f0afL, 0xeee2ed18L, 0xf0a5bd1dL, 0xf464a0aaL, 0xf9278673L, 0xfde69bc4L, 0x89b8fd09L, 0x8d79e0beL, 0x803ac667L, 0x84fbdbd0L, 0x9abc8bd5L, 0x9e7d9662L, 0x933eb0bbL, 0x97ffad0cL, 0xafb010b1L, 0xab710d06L, 0xa6322bdfL, 0xa2f33668L, 0xbcb4666dL, 0xb8757bdaL, 0xb5365d03L, 0xb1f740b4L};
/*-- Ugly, innit? --*/

/*-------------------------------------------------------------*/

/*--- end                                        crctable.c ---*/

/*-------------------------------------------------------------*/

/*-------------------------------------------------------------*/

/*--- Table for randomising repetitive blocks               ---*/

/*---                                           randtable.c ---*/

/*-------------------------------------------------------------*/

/*---------------------------------------------*/

Int32 BZ2_rNums[512] = {619, 720, 127, 481, 931, 816, 813, 233, 566, 247, 985, 724, 205, 454, 863, 491, 741, 242, 949, 214, 733, 859, 335, 708, 621, 574, 73, 654, 730, 472, 419, 436, 278, 496, 867, 210, 399, 680, 480, 51, 878, 465, 811, 169, 869, 675, 611, 697, 867, 561, 862, 687, 507, 283, 482, 129, 807, 591, 733, 623, 150, 238, 59, 379, 684, 877, 625, 169, 643, 105, 170, 607, 520, 932, 727, 476, 693, 425, 174, 647, 73, 122, 335, 530, 442, 853, 695, 249, 445, 515, 909, 545, 703, 919, 874, 474, 882, 500, 594, 612, 641, 801, 220, 162, 819, 984, 589, 513, 495, 799, 161, 604, 958, 533, 221, 400, 386, 867, 600, 782, 382, 596, 414, 171, 516, 375, 682, 485, 911, 276, 98, 553, 163, 354, 666, 933, 424, 341, 533, 870, 227, 730, 475, 186, 263, 647, 537, 686, 600, 224, 469, 68, 770, 919, 190, 373, 294, 822, 808, 206, 184, 943, 795, 384, 383, 461, 404, 758, 839, 887, 715, 67, 618, 276, 204, 918, 873, 777, 604, 560, 951, 160, 578, 722, 79, 804, 96, 409, 713, 940, 652, 934, 970, 447, 318, 353, 859, 672, 112, 785, 645, 863, 803, 350, 139, 93, 354, 99, 820, 908, 609, 772, 154, 274, 580, 184, 79, 626, 630, 742, 653, 282, 762, 623, 680, 81, 927, 626, 789, 125, 411, 521, 938, 300, 821, 78, 343, 175, 128, 250, 170, 774, 972, 275, 999, 639, 495, 78, 352, 126, 857, 956, 358, 619, 580, 124, 737, 594, 701, 612, 669, 112, 134, 694, 363, 992, 809, 743, 168, 974, 944, 375, 748, 52, 600, 747, 642, 182, 862, 81, 344, 805, 988, 739, 511, 655, 814, 334, 249, 515, 897, 955, 664, 981, 649, 113, 974, 459, 893, 228, 433, 837, 553, 268, 926, 240, 102, 654, 459, 51, 686, 754, 806, 760, 493, 403, 415, 394, 687, 700, 946, 670, 656, 610, 738, 392, 760, 799, 887, 653, 978, 321, 576, 617, 626, 502, 894, 679, 243, 440, 680, 879, 194, 572, 640, 724, 926, 56, 204, 700, 707, 151, 457, 449, 797, 195, 791, 558, 945, 679, 297, 59, 87, 824, 713, 663, 412, 693, 342, 606, 134, 108, 571, 364, 631, 212, 174, 643, 304, 329, 343, 97, 430, 751, 497, 314, 983, 374, 822, 928, 140, 206, 73, 263, 980, 736, 876, 478, 430, 305, 170, 514, 364, 692, 829, 82, 855, 953, 676, 246, 369, 970, 294, 750, 807, 827, 150, 790, 288, 923, 804, 378, 215, 828, 592, 281, 565, 555, 710, 82, 896, 831, 547, 261, 524, 462, 293, 465, 502, 56, 661, 821, 976, 991, 658, 869, 905, 758, 745, 193, 768, 550, 608, 933, 378, 286, 215, 979, 792, 961, 61, 688, 793, 644, 986, 403, 106, 366, 905, 644, 372, 567, 466, 434, 645, 210, 389, 550, 919, 135, 780, 773, 635, 389, 707, 100, 626, 958, 165, 504, 920, 176, 193, 713, 857, 265, 203, 50, 668, 108, 645, 990, 626, 197, 510, 357, 358, 850, 858, 364, 936, 638};
/*-------------------------------------------------------------*/
/*--- end                                       randtable.c ---*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*--- Compression machinery (not incl block sorting)        ---*/
/*---                                            compress.c ---*/
/*-------------------------------------------------------------*/
/*---------------------------------------------------*/
/*--- Bit stream I/O                              ---*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
void BZ2_bsInitWrite(EState *s) {
   s->bsLive = 0;
   s->bsBuff = 0;
}

/*---------------------------------------------------*/
static void bsFinishWrite(EState *s) {
   while(s->bsLive > 0) {
      s->zbits[s->numZ] = (UChar) (s->bsBuff >> 24);
      s->numZ++;
      s->bsBuff <<= 8;
      s->bsLive -= 8;
   }
}

/*---------------------------------------------------*/
/*---------------------------------------------------*/
inline static void bsW(EState *s, Int32 n, UInt32 v) {
   {
      while(s->bsLive >= 8) {
         s->zbits[s->numZ] = (UChar) (s->bsBuff >> 24);
         s->numZ++;
         s->bsBuff <<= 8;
         s->bsLive -= 8;
      }
   }
   ;
   s->bsBuff |= (v << (32 - s->bsLive - n));
   s->bsLive += n;
}

/*---------------------------------------------------*/
static void bsPutUInt32(EState *s, UInt32 u) {
   clava_dcg_global[ 17 ]++;
   bsW(s, 8, (u >> 24) & 0xffL);
   clava_dcg_global[ 17 ]++;
   bsW(s, 8, (u >> 16) & 0xffL);
   clava_dcg_global[ 17 ]++;
   bsW(s, 8, (u >> 8) & 0xffL);
   clava_dcg_global[ 17 ]++;
   bsW(s, 8, u & 0xffL);
}

/*---------------------------------------------------*/
static void bsPutUChar(EState *s, UChar c) {
   clava_dcg_global[ 18 ]++;
   bsW(s, 8, (UInt32) c);
}

/*---------------------------------------------------*/
/*--- The back end proper                         ---*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
static void makeMaps_e(EState *s) {
   Int32 i;
   s->nInUse = 0;
   for(i = 0; i < 256; i++) if(s->inUse[i]) {
      s->unseqToSeq[i] = s->nInUse;
      s->nInUse++;
   }
}

/*---------------------------------------------------*/
static void generateMTFValues(EState *s) {
   UChar yy[256];
   Int32 i, j;
   Int32 zPend;
   Int32 wr;
   Int32 EOB;
   /*
   After sorting (eg, here),
   s->arr1 [ 0 .. s->nblock-1 ] holds sorted order,
   and
   ((UChar*)s->arr2) [ 0 .. s->nblock-1 ]
   holds the original block data.
   
   The first thing to do is generate the MTF values,
   and put them in
   ((UInt16*)s->arr1) [ 0 .. s->nblock-1 ].
   Because there are strictly fewer or equal MTF values
   than block values, ptr values in this area are overwritten
   with MTF values only when they are no longer needed.
   
   The final compressed bitstream is generated into the
   area starting at
   (UChar*) (&((UChar*)s->arr2)[s->nblock])
   
   These storage aliases are set up in bzCompressInit(),
   except for the last one, which is arranged in
   compressBlock().
   */
   UInt32 *ptr = s->ptr;
   UChar *block = s->block;
   UInt16 *mtfv = s->mtfv;
   clava_dcg_global[ 19 ]++;
   makeMaps_e(s);
   EOB = s->nInUse + 1;
   for(i = 0; i <= EOB; i++) s->mtfFreq[i] = 0;
   wr = 0;
   zPend = 0;
   for(i = 0; i < s->nInUse; i++) yy[i] = (UChar) i;
   for(i = 0; i < s->nblock; i++) {
      UChar ll_i;
      ;
      j = ptr[i] - 1;
      if(j < 0) j += s->nblock;
      ll_i = s->unseqToSeq[block[j]];
      ;
      if(yy[0] == ll_i) {
         zPend++;
      }
      else {
         if(zPend > 0) {
            zPend--;
            while(((Bool) 1)) {
               if(zPend & 1) {
                  mtfv[wr] = 1;
                  wr++;
                  s->mtfFreq[1]++;
               }
               else {
                  mtfv[wr] = 0;
                  wr++;
                  s->mtfFreq[0]++;
               }
               if(zPend < 2) break;
               zPend = (zPend - 2) / 2;
            }
            ;
            zPend = 0;
         }
         {
            register UChar rtmp;
            register UChar *ryy_j;
            register UChar rll_i;
            rtmp = yy[1];
            yy[1] = yy[0];
            ryy_j = &(yy[1]);
            rll_i = ll_i;
            while(rll_i != rtmp) {
               register UChar rtmp2;
               ryy_j++;
               rtmp2 = rtmp;
               rtmp = *ryy_j;
               *ryy_j = rtmp2;
            }
            ;
            yy[0] = rtmp;
            j = ryy_j - &(yy[0]);
            mtfv[wr] = j + 1;
            wr++;
            s->mtfFreq[j + 1]++;
         }
      }
   }
   if(zPend > 0) {
      zPend--;
      while(((Bool) 1)) {
         if(zPend & 1) {
            mtfv[wr] = 1;
            wr++;
            s->mtfFreq[1]++;
         }
         else {
            mtfv[wr] = 0;
            wr++;
            s->mtfFreq[0]++;
         }
         if(zPend < 2) break;
         zPend = (zPend - 2) / 2;
      }
      ;
      zPend = 0;
   }
   mtfv[wr] = EOB;
   wr++;
   s->mtfFreq[EOB]++;
   s->nMTF = wr;
}

/*---------------------------------------------------*/
static void sendMTFValues(EState *s) {
   Int32 v, t, i, j, gs, ge, totc, bt, bc, iter;
   Int32 nSelectors, alphaSize, minLen, maxLen, selCtr;
   Int32 nGroups, nBytes;
   /*--
   UChar  len [BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
   is a global since the decoder also needs it.
   
   Int32  code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
   Int32  rfreq[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
   are also globals only used in this proc.
   Made global to keep stack frame size small.
   --*/
   UInt16 cost[6];
   Int32 fave[6];
   UInt16 *mtfv = s->mtfv;
   if(s->verbosity >= 3) {
      clava_dcg_global[ 20 ]++;
      fprintf(stderr, "      %d in block, %d after MTF & 1-2 coding, %d+2 syms in use\n", s->nblock, s->nMTF, s->nInUse);
   }
   alphaSize = s->nInUse + 2;
   for(t = 0; t < 6; t++) for(v = 0; v < alphaSize; v++) s->len[t][v] = 15;
   /*--- Decide how many coding tables to use ---*/
   {
      if(!(s->nMTF > 0)) {
         clava_dcg_global[ 21 ]++;
         BZ2_bz__AssertH__fail(3001);
      }
   }
   ;
   if(s->nMTF < 200) nGroups = 2;
   else if(s->nMTF < 600) nGroups = 3;
   else if(s->nMTF < 1200) nGroups = 4;
   else if(s->nMTF < 2400) nGroups = 5;
   else nGroups = 6;
   /*--- Generate an initial set of coding tables ---*/
   {
      Int32 nPart, remF, tFreq, aFreq;
      nPart = nGroups;
      remF = s->nMTF;
      gs = 0;
      while(nPart > 0) {
         tFreq = remF / nPart;
         ge = gs - 1;
         aFreq = 0;
         while(aFreq < tFreq && ge < alphaSize - 1) {
            ge++;
            aFreq += s->mtfFreq[ge];
         }
         if(ge > gs && nPart != nGroups && nPart != 1 && ((nGroups - nPart) % 2 == 1)) {
            aFreq -= s->mtfFreq[ge];
            ge--;
         }
         if(s->verbosity >= 3) {
            clava_dcg_global[ 20 ]++;
            fprintf(stderr, "      initial group %d, [%d .. %d], has %d syms (%4.1f%%)\n", nPart, gs, ge, aFreq, (100.0 * (float) aFreq) / (float) (s->nMTF));
         }
         for(v = 0; v < alphaSize; v++) if(v >= gs && v <= ge) s->len[nPart - 1][v] = 0;
         else s->len[nPart - 1][v] = 15;
         nPart--;
         gs = ge + 1;
         remF -= aFreq;
      }
   }
   /*---
   Iterate up to BZ_N_ITERS times to improve the tables.
   ---*/
   for(iter = 0; iter < 4; iter++) {
      for(t = 0; t < nGroups; t++) fave[t] = 0;
      for(t = 0; t < nGroups; t++) for(v = 0; v < alphaSize; v++) s->rfreq[t][v] = 0;
      /*---
      Set up an auxiliary length table which is used to fast-track
      the common case (nGroups == 6).
      ---*/
      if(nGroups == 6) {
         for(v = 0; v < alphaSize; v++) {
            s->len_pack[v][0] = (s->len[1][v] << 16) | s->len[0][v];
            s->len_pack[v][1] = (s->len[3][v] << 16) | s->len[2][v];
            s->len_pack[v][2] = (s->len[5][v] << 16) | s->len[4][v];
         }
      }
      nSelectors = 0;
      totc = 0;
      gs = 0;
      while(((Bool) 1)) {
         /*--- Set group start & end marks. --*/
         if(gs >= s->nMTF) break;
         ge = gs + 50 - 1;
         if(ge >= s->nMTF) ge = s->nMTF - 1;
         /*--
         Calculate the cost of this group as coded
         by each of the coding tables.
         --*/
         for(t = 0; t < nGroups; t++) cost[t] = 0;
         if(nGroups == 6 && 50 == ge - gs + 1) {
            /*--- fast track the common case ---*/
            register UInt32 cost01, cost23, cost45;
            register UInt16 icv;
            cost01 = cost23 = cost45 = 0;
            icv = mtfv[gs + (0)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (1)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (2)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (3)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (4)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (5)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (6)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (7)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (8)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (9)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (10)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (11)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (12)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (13)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (14)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (15)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (16)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (17)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (18)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (19)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (20)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (21)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (22)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (23)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (24)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (25)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (26)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (27)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (28)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (29)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (30)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (31)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (32)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (33)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (34)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (35)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (36)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (37)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (38)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (39)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (40)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (41)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (42)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (43)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (44)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (45)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (46)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (47)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (48)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            icv = mtfv[gs + (49)];
            cost01 += s->len_pack[icv][0];
            cost23 += s->len_pack[icv][1];
            cost45 += s->len_pack[icv][2];
            ;
            cost[0] = cost01 & 0xffff;
            cost[1] = cost01 >> 16;
            cost[2] = cost23 & 0xffff;
            cost[3] = cost23 >> 16;
            cost[4] = cost45 & 0xffff;
            cost[5] = cost45 >> 16;
         }
         else {
            /*--- slow version which correctly handles all situations ---*/
            for(i = gs; i <= ge; i++) {
               UInt16 icv = mtfv[i];
               for(t = 0; t < nGroups; t++) cost[t] += s->len[t][icv];
            }
         }
         /*--
         Find the coding table which is best for this group,
         and record its identity in the selector table.
         --*/
         bc = 999999999;
         bt = -1;
         for(t = 0; t < nGroups; t++) if(cost[t] < bc) {
            bc = cost[t];
            bt = t;
         }
         ;
         totc += bc;
         fave[bt]++;
         s->selector[nSelectors] = bt;
         nSelectors++;
         /*--
         Increment the symbol frequencies for the selected table.
         --*/
         if(nGroups == 6 && 50 == ge - gs + 1) {
            /*--- fast track the common case ---*/
            s->rfreq[bt][mtfv[gs + (0)]]++;
            s->rfreq[bt][mtfv[gs + (1)]]++;
            s->rfreq[bt][mtfv[gs + (2)]]++;
            s->rfreq[bt][mtfv[gs + (3)]]++;
            s->rfreq[bt][mtfv[gs + (4)]]++;
            s->rfreq[bt][mtfv[gs + (5)]]++;
            s->rfreq[bt][mtfv[gs + (6)]]++;
            s->rfreq[bt][mtfv[gs + (7)]]++;
            s->rfreq[bt][mtfv[gs + (8)]]++;
            s->rfreq[bt][mtfv[gs + (9)]]++;
            s->rfreq[bt][mtfv[gs + (10)]]++;
            s->rfreq[bt][mtfv[gs + (11)]]++;
            s->rfreq[bt][mtfv[gs + (12)]]++;
            s->rfreq[bt][mtfv[gs + (13)]]++;
            s->rfreq[bt][mtfv[gs + (14)]]++;
            s->rfreq[bt][mtfv[gs + (15)]]++;
            s->rfreq[bt][mtfv[gs + (16)]]++;
            s->rfreq[bt][mtfv[gs + (17)]]++;
            s->rfreq[bt][mtfv[gs + (18)]]++;
            s->rfreq[bt][mtfv[gs + (19)]]++;
            s->rfreq[bt][mtfv[gs + (20)]]++;
            s->rfreq[bt][mtfv[gs + (21)]]++;
            s->rfreq[bt][mtfv[gs + (22)]]++;
            s->rfreq[bt][mtfv[gs + (23)]]++;
            s->rfreq[bt][mtfv[gs + (24)]]++;
            s->rfreq[bt][mtfv[gs + (25)]]++;
            s->rfreq[bt][mtfv[gs + (26)]]++;
            s->rfreq[bt][mtfv[gs + (27)]]++;
            s->rfreq[bt][mtfv[gs + (28)]]++;
            s->rfreq[bt][mtfv[gs + (29)]]++;
            s->rfreq[bt][mtfv[gs + (30)]]++;
            s->rfreq[bt][mtfv[gs + (31)]]++;
            s->rfreq[bt][mtfv[gs + (32)]]++;
            s->rfreq[bt][mtfv[gs + (33)]]++;
            s->rfreq[bt][mtfv[gs + (34)]]++;
            s->rfreq[bt][mtfv[gs + (35)]]++;
            s->rfreq[bt][mtfv[gs + (36)]]++;
            s->rfreq[bt][mtfv[gs + (37)]]++;
            s->rfreq[bt][mtfv[gs + (38)]]++;
            s->rfreq[bt][mtfv[gs + (39)]]++;
            s->rfreq[bt][mtfv[gs + (40)]]++;
            s->rfreq[bt][mtfv[gs + (41)]]++;
            s->rfreq[bt][mtfv[gs + (42)]]++;
            s->rfreq[bt][mtfv[gs + (43)]]++;
            s->rfreq[bt][mtfv[gs + (44)]]++;
            s->rfreq[bt][mtfv[gs + (45)]]++;
            s->rfreq[bt][mtfv[gs + (46)]]++;
            s->rfreq[bt][mtfv[gs + (47)]]++;
            s->rfreq[bt][mtfv[gs + (48)]]++;
            s->rfreq[bt][mtfv[gs + (49)]]++;
         }
         else {
            /*--- slow version which correctly handles all situations ---*/
            for(i = gs; i <= ge; i++) s->rfreq[bt][mtfv[i]]++;
         }
         gs = ge + 1;
      }
      if(s->verbosity >= 3) {
         clava_dcg_global[ 20 ]++;
         fprintf(stderr, "      pass %d: size is %d, grp uses are ", iter + 1, totc / 8);
         for(t = 0; t < nGroups; t++) {
            clava_dcg_global[ 20 ]++;
            fprintf(stderr, "%d ", fave[t]);
         }
         clava_dcg_global[ 20 ]++;
         fprintf(stderr, "\n");
      }
      /*--
      Recompute the tables based on the accumulated frequencies.
      --*/
      for(t = 0; t < nGroups; t++) {
         clava_dcg_global[ 22 ]++;
         BZ2_hbMakeCodeLengths(&(s->len[t][0]), &(s->rfreq[t][0]), alphaSize, 20);
      }
   }
   {
      if(!(nGroups < 8)) {
         clava_dcg_global[ 21 ]++;
         BZ2_bz__AssertH__fail(3002);
      }
   }
   ;
   {
      if(!(nSelectors < 32768 && nSelectors <= (2 + (900000 / 50)))) {
         clava_dcg_global[ 21 ]++;
         BZ2_bz__AssertH__fail(3003);
      }
   }
   ;
   /*--- Compute MTF values for the selectors. ---*/
   {
      UChar pos[6];
      UChar ll_i;
      UChar tmp2;
      UChar tmp;
      for(i = 0; i < nGroups; i++) pos[i] = i;
      for(i = 0; i < nSelectors; i++) {
         ll_i = s->selector[i];
         j = 0;
         tmp = pos[j];
         while(ll_i != tmp) {
            j++;
            tmp2 = tmp;
            tmp = pos[j];
            pos[j] = tmp2;
         }
         ;
         pos[0] = tmp;
         s->selectorMtf[i] = j;
      }
   }
   ;
   /*--- Assign actual codes for the tables. --*/
   for(t = 0; t < nGroups; t++) {
      minLen = 32;
      maxLen = 0;
      for(i = 0; i < alphaSize; i++) {
         if(s->len[t][i] > maxLen) maxLen = s->len[t][i];
         if(s->len[t][i] < minLen) minLen = s->len[t][i];
      }
      {
         if(!(!(maxLen > 20))) {
            clava_dcg_global[ 21 ]++;
            BZ2_bz__AssertH__fail(3004);
         }
      }
      ;
      {
         if(!(!(minLen < 1))) {
            clava_dcg_global[ 21 ]++;
            BZ2_bz__AssertH__fail(3005);
         }
      }
      ;
      clava_dcg_global[ 23 ]++;
      BZ2_hbAssignCodes(&(s->code[t][0]), &(s->len[t][0]), minLen, maxLen, alphaSize);
   }
   /*--- Transmit the mapping table. ---*/
   {
      Bool inUse16[16];
      for(i = 0; i < 16; i++) {
         inUse16[i] = ((Bool) 0);
         for(j = 0; j < 16; j++) if(s->inUse[i * 16 + j]) inUse16[i] = ((Bool) 1);
      }
      nBytes = s->numZ;
      for(i = 0; i < 16; i++) if(inUse16[i]) {
         clava_dcg_global[ 24 ]++;
         bsW(s, 1, 1);
      }
      else {
         clava_dcg_global[ 24 ]++;
         bsW(s, 1, 0);
      }
      for(i = 0; i < 16; i++) if(inUse16[i]) for(j = 0; j < 16; j++) {
         if(s->inUse[i * 16 + j]) {
            clava_dcg_global[ 24 ]++;
            bsW(s, 1, 1);
         }
         else {
            clava_dcg_global[ 24 ]++;
            bsW(s, 1, 0);
         }
      }
      if(s->verbosity >= 3) {
         clava_dcg_global[ 20 ]++;
         fprintf(stderr, "      bytes: mapping %d, ", s->numZ - nBytes);
      }
   }
   /*--- Now the selectors. ---*/
   nBytes = s->numZ;
   clava_dcg_global[ 24 ]++;
   bsW(s, 3, nGroups);
   clava_dcg_global[ 24 ]++;
   bsW(s, 15, nSelectors);
   for(i = 0; i < nSelectors; i++) {
      for(j = 0; j < s->selectorMtf[i]; j++) {
         clava_dcg_global[ 24 ]++;
         bsW(s, 1, 1);
      }
      clava_dcg_global[ 24 ]++;
      bsW(s, 1, 0);
   }
   if(s->verbosity >= 3) {
      clava_dcg_global[ 20 ]++;
      fprintf(stderr, "selectors %d, ", s->numZ - nBytes);
   }
   /*--- Now the coding tables. ---*/
   nBytes = s->numZ;
   for(t = 0; t < nGroups; t++) {
      Int32 curr = s->len[t][0];
      clava_dcg_global[ 24 ]++;
      bsW(s, 5, curr);
      for(i = 0; i < alphaSize; i++) {
         while(curr < s->len[t][i]) {
            clava_dcg_global[ 24 ]++;
            bsW(s, 2, 2);
            curr++;
         }
         ;
         /*10*/
         while(curr > s->len[t][i]) {
            clava_dcg_global[ 24 ]++;
            bsW(s, 2, 3);
            curr--;
         }
         ;
         /*11*/
         clava_dcg_global[ 24 ]++;
         bsW(s, 1, 0);
      }
   }
   if(s->verbosity >= 3) {
      clava_dcg_global[ 20 ]++;
      fprintf(stderr, "code lengths %d, ", s->numZ - nBytes);
   }
   /*--- And finally, the block data proper ---*/
   nBytes = s->numZ;
   selCtr = 0;
   gs = 0;
   while(((Bool) 1)) {
      if(gs >= s->nMTF) break;
      ge = gs + 50 - 1;
      if(ge >= s->nMTF) ge = s->nMTF - 1;
      {
         if(!(s->selector[selCtr] < nGroups)) {
            clava_dcg_global[ 21 ]++;
            BZ2_bz__AssertH__fail(3006);
         }
      }
      ;
      if(nGroups == 6 && 50 == ge - gs + 1) {
         /*--- fast track the common case ---*/
         UInt16 mtfv_i;
         UChar *s_len_sel_selCtr = &(s->len[s->selector[selCtr]][0]);
         Int32 *s_code_sel_selCtr = &(s->code[s->selector[selCtr]][0]);
         mtfv_i = mtfv[gs + (0)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (1)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (2)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (3)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (4)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (5)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (6)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (7)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (8)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (9)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (10)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (11)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (12)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (13)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (14)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (15)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (16)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (17)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (18)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (19)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (20)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (21)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (22)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (23)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (24)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (25)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (26)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (27)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (28)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (29)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (30)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (31)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (32)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (33)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (34)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (35)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (36)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (37)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (38)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (39)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (40)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (41)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (42)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (43)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (44)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (45)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (46)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (47)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (48)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
         mtfv_i = mtfv[gs + (49)];
         clava_dcg_global[ 24 ]++;
         bsW(s, s_len_sel_selCtr[mtfv_i], s_code_sel_selCtr[mtfv_i]);
      }
      else {
         /*--- slow version which correctly handles all situations ---*/
         for(i = gs; i <= ge; i++) {
            clava_dcg_global[ 24 ]++;
            bsW(s, s->len[s->selector[selCtr]][mtfv[i]], s->code[s->selector[selCtr]][mtfv[i]]);
         }
      }
      gs = ge + 1;
      selCtr++;
   }
   {
      if(!(selCtr == nSelectors)) {
         clava_dcg_global[ 21 ]++;
         BZ2_bz__AssertH__fail(3007);
      }
   }
   ;
   if(s->verbosity >= 3) {
      clava_dcg_global[ 20 ]++;
      fprintf(stderr, "codes %d\n", s->numZ - nBytes);
   }
}

/*---------------------------------------------------*/
void BZ2_compressBlock(EState *s, Bool is_last_block) {
   if(s->nblock > 0) {
      {
         s->blockCRC = ~(s->blockCRC);
      }
      ;
      s->combinedCRC = (s->combinedCRC << 1) | (s->combinedCRC >> 31);
      s->combinedCRC ^= s->blockCRC;
      if(s->blockNo > 1) s->numZ = 0;
      if(s->verbosity >= 2) {
         clava_dcg_global[ 25 ]++;
         fprintf(stderr, "    block %d: crc = 0x%8x, combined CRC = 0x%8x, size = %d\n", s->blockNo, s->blockCRC, s->combinedCRC, s->nblock);
      }
      clava_dcg_global[ 26 ]++;
      BZ2_blockSort(s);
   }
   s->zbits = (UChar *) (&((UChar *) s->arr2)[s->nblock]);
   /*-- If this is the first block, create the stream header. --*/
   if(s->blockNo == 1) {
      clava_dcg_global[ 27 ]++;
      BZ2_bsInitWrite(s);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x42);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x5a);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x68);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, (UChar) (0x30 + s->blockSize100k));
   }
   if(s->nblock > 0) {
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x31);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x41);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x59);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x26);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x53);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x59);
      /*-- Now the block's CRC, so it is in a known place. --*/
      clava_dcg_global[ 29 ]++;
      bsPutUInt32(s, s->blockCRC);
      /*--
      Now a single bit indicating (non-)randomisation.
      As of version 0.9.5, we use a better sorting algorithm
      which makes randomisation unnecessary.  So always set
      the randomised bit to 'no'.  Of course, the decoder
      still needs to be able to handle randomised blocks
      so as to maintain backwards compatibility with
      older versions of bzip2.
      --*/
      clava_dcg_global[ 30 ]++;
      bsW(s, 1, 0);
      clava_dcg_global[ 30 ]++;
      bsW(s, 24, s->origPtr);
      clava_dcg_global[ 31 ]++;
      generateMTFValues(s);
      clava_dcg_global[ 32 ]++;
      sendMTFValues(s);
   }
   /*-- If this is the last block, add the stream trailer. --*/
   if(is_last_block) {
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x17);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x72);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x45);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x38);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x50);
      clava_dcg_global[ 28 ]++;
      bsPutUChar(s, 0x90);
      clava_dcg_global[ 29 ]++;
      bsPutUInt32(s, s->combinedCRC);
      if(s->verbosity >= 2) {
         clava_dcg_global[ 25 ]++;
         fprintf(stderr, "    final combined CRC = 0x%x\n   ", s->combinedCRC);
      }
      clava_dcg_global[ 33 ]++;
      bsFinishWrite(s);
   }
}

/*-------------------------------------------------------------*/
/*--- end                                        compress.c ---*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*--- Decompression machinery                               ---*/
/*---                                          decompress.c ---*/
/*-------------------------------------------------------------*/
/*---------------------------------------------------*/
static void makeMaps_d(DState *s) {
   Int32 i;
   s->nInUse = 0;
   for(i = 0; i < 256; i++) if(s->inUse[i]) {
      s->seqToUnseq[s->nInUse] = i;
      s->nInUse++;
   }
}

/*---------------------------------------------------*/
/*---------------------------------------------------*/
/*the longest code*/
/*---------------------------------------------------*/
Int32 BZ2_decompress(DState *s) {
   UChar uc;
   Int32 retVal;
   Int32 minLen, maxLen;
   bz_stream *strm = s->strm;
   /*stuff that needs to be saved/restored*/
   Int32 i;
   Int32 j;
   Int32 t;
   Int32 alphaSize;
   Int32 nGroups;
   Int32 nSelectors;
   Int32 EOB;
   Int32 groupNo;
   Int32 groupPos;
   Int32 nextSym;
   Int32 nblockMAX;
   Int32 nblock;
   Int32 es;
   Int32 N;
   Int32 curr;
   Int32 zt;
   Int32 zn;
   Int32 zvec;
   Int32 zj;
   Int32 gSel;
   Int32 gMinlen;
   Int32 *gLimit;
   Int32 *gBase;
   Int32 *gPerm;
   if(s->state == 10) {
      /*initialise the save area*/
      s->save_i = 0;
      s->save_j = 0;
      s->save_t = 0;
      s->save_alphaSize = 0;
      s->save_nGroups = 0;
      s->save_nSelectors = 0;
      s->save_EOB = 0;
      s->save_groupNo = 0;
      s->save_groupPos = 0;
      s->save_nextSym = 0;
      s->save_nblockMAX = 0;
      s->save_nblock = 0;
      s->save_es = 0;
      s->save_N = 0;
      s->save_curr = 0;
      s->save_zt = 0;
      s->save_zn = 0;
      s->save_zvec = 0;
      s->save_zj = 0;
      s->save_gSel = 0;
      s->save_gMinlen = 0;
      s->save_gLimit = ((void *) 0);
      s->save_gBase = ((void *) 0);
      s->save_gPerm = ((void *) 0);
   }
   /*restore from the save area*/
   i = s->save_i;
   j = s->save_j;
   t = s->save_t;
   alphaSize = s->save_alphaSize;
   nGroups = s->save_nGroups;
   nSelectors = s->save_nSelectors;
   EOB = s->save_EOB;
   groupNo = s->save_groupNo;
   groupPos = s->save_groupPos;
   nextSym = s->save_nextSym;
   nblockMAX = s->save_nblockMAX;
   nblock = s->save_nblock;
   es = s->save_es;
   N = s->save_N;
   curr = s->save_curr;
   zt = s->save_zt;
   zn = s->save_zn;
   zvec = s->save_zvec;
   zj = s->save_zj;
   gSel = s->save_gSel;
   gMinlen = s->save_gMinlen;
   gLimit = s->save_gLimit;
   gBase = s->save_gBase;
   gPerm = s->save_gPerm;
   retVal = 0;
   switch (s->state) {
      case 10:
      s->state = 10;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x42) {
         retVal = (-5);
         goto save_state_and_return;
      }
      ;
      ;
      case 11:
      s->state = 11;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x5a) {
         retVal = (-5);
         goto save_state_and_return;
      }
      ;
      ;
      case 12:
      s->state = 12;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      if(uc != 0x68) {
         retVal = (-5);
         goto save_state_and_return;
      }
      ;
      ;
      case 13:
      s->state = 13;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            s->blockSize100k = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      if(s->blockSize100k < (0x30 + 1) || s->blockSize100k > (0x30 + 9)) {
         retVal = (-5);
         goto save_state_and_return;
      }
      ;
      ;
      s->blockSize100k -= 0x30;
      if(s->smallDecompress) {
         clava_dcg_global[ 34 ]++;
         s->ll16 = (strm->bzalloc)(strm->opaque, (s->blockSize100k * 100000 * sizeof(UInt16)), 1);
         clava_dcg_global[ 34 ]++;
         s->ll4 = (strm->bzalloc)(strm->opaque, (((1 + s->blockSize100k * 100000) >> 1) * sizeof(UChar)), 1);
         if(s->ll16 == ((void *) 0) || s->ll4 == ((void *) 0)) {
            retVal = (-3);
            goto save_state_and_return;
         }
         ;
         ;
      }
      else {
         clava_dcg_global[ 34 ]++;
         s->tt = (strm->bzalloc)(strm->opaque, (s->blockSize100k * 100000 * sizeof(Int32)), 1);
         if(s->tt == ((void *) 0)) {
            retVal = (-3);
            goto save_state_and_return;
         }
         ;
         ;
      }
      case 14:
      s->state = 14;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc == 0x17) goto endhdr_2;
      if(uc != 0x31) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 15:
      s->state = 15;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x41) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 16:
      s->state = 16;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x59) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 17:
      s->state = 17;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x26) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 18:
      s->state = 18;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x53) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 19:
      s->state = 19;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x59) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      s->currBlockNo++;
      if(s->verbosity >= 2) {
         clava_dcg_global[ 35 ]++;
         fprintf(stderr, "\n    [%d: huff+mtf ", s->currBlockNo);
      }
      s->storedBlockCRC = 0;
      case 20:
      s->state = 20;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32) uc);
      case 21:
      s->state = 21;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32) uc);
      case 22:
      s->state = 22;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32) uc);
      case 23:
      s->state = 23;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32) uc);
      case 24:
      s->state = 24;
      while(((Bool) 1)) {
         if(s->bsLive >= 1) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
            s->bsLive -= 1;
            s->blockRandomised = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->origPtr = 0;
      case 25:
      s->state = 25;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->origPtr = (s->origPtr << 8) | ((Int32) uc);
      case 26:
      s->state = 26;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->origPtr = (s->origPtr << 8) | ((Int32) uc);
      case 27:
      s->state = 27;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->origPtr = (s->origPtr << 8) | ((Int32) uc);
      if(s->origPtr < 0) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      if(s->origPtr > 10 + 100000 * s->blockSize100k) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      /*--- Receive the mapping table ---*/
      for(i = 0; i < 16; i++) {
         case 28:
         s->state = 28;
         while(((Bool) 1)) {
            if(s->bsLive >= 1) {
               UInt32 v;
               v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
               s->bsLive -= 1;
               uc = v;
               break;
            }
            if(s->strm->avail_in == 0) {
               retVal = 0;
               goto save_state_and_return;
            }
            ;
            ;
            s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
            s->bsLive += 8;
            s->strm->next_in++;
            s->strm->avail_in--;
            s->strm->total_in_lo32++;
            if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
         }
         ;
         if(uc == 1) s->inUse16[i] = ((Bool) 1);
         else s->inUse16[i] = ((Bool) 0);
      }
      for(i = 0; i < 256; i++) s->inUse[i] = ((Bool) 0);
      for(i = 0; i < 16; i++) if(s->inUse16[i]) for(j = 0; j < 16; j++) {
         case 29:
         s->state = 29;
         while(((Bool) 1)) {
            if(s->bsLive >= 1) {
               UInt32 v;
               v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
               s->bsLive -= 1;
               uc = v;
               break;
            }
            if(s->strm->avail_in == 0) {
               retVal = 0;
               goto save_state_and_return;
            }
            ;
            ;
            s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
            s->bsLive += 8;
            s->strm->next_in++;
            s->strm->avail_in--;
            s->strm->total_in_lo32++;
            if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
         }
         ;
         if(uc == 1) s->inUse[i * 16 + j] = ((Bool) 1);
      }
      clava_dcg_global[ 36 ]++;
      makeMaps_d(s);
      if(s->nInUse == 0) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      alphaSize = s->nInUse + 2;
      /*--- Now the selectors ---*/
      case 30:
      s->state = 30;
      while(((Bool) 1)) {
         if(s->bsLive >= 3) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 3)) & ((1 << 3) - 1);
            s->bsLive -= 3;
            nGroups = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(nGroups < 2 || nGroups > 6) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 31:
      s->state = 31;
      while(((Bool) 1)) {
         if(s->bsLive >= 15) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 15)) & ((1 << 15) - 1);
            s->bsLive -= 15;
            nSelectors = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(nSelectors < 1) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      for(i = 0; i < nSelectors; i++) {
         j = 0;
         while(((Bool) 1)) {
            case 32:
            s->state = 32;
            while(((Bool) 1)) {
               if(s->bsLive >= 1) {
                  UInt32 v;
                  v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                  s->bsLive -= 1;
                  uc = v;
                  break;
               }
               if(s->strm->avail_in == 0) {
                  retVal = 0;
                  goto save_state_and_return;
               }
               ;
               ;
               s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
               s->bsLive += 8;
               s->strm->next_in++;
               s->strm->avail_in--;
               s->strm->total_in_lo32++;
               if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
            }
            ;
            if(uc == 0) break;
            j++;
            if(j >= nGroups) {
               retVal = (-4);
               goto save_state_and_return;
            }
            ;
            ;
         }
         s->selectorMtf[i] = j;
      }
      /*--- Undo the MTF values for the selectors. ---*/
      {
         UChar pos[6];
         UChar tmp;
         UChar v;
         for(v = 0; v < nGroups; v++) pos[v] = v;
         for(i = 0; i < nSelectors; i++) {
            v = s->selectorMtf[i];
            tmp = pos[v];
            while(v > 0) {
               pos[v] = pos[v - 1];
               v--;
            }
            pos[0] = tmp;
            s->selector[i] = tmp;
         }
      }
      /*--- Now the coding tables ---*/
      for(t = 0; t < nGroups; t++) {
         case 33:
         s->state = 33;
         while(((Bool) 1)) {
            if(s->bsLive >= 5) {
               UInt32 v;
               v = (s->bsBuff >> (s->bsLive - 5)) & ((1 << 5) - 1);
               s->bsLive -= 5;
               curr = v;
               break;
            }
            if(s->strm->avail_in == 0) {
               retVal = 0;
               goto save_state_and_return;
            }
            ;
            ;
            s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
            s->bsLive += 8;
            s->strm->next_in++;
            s->strm->avail_in--;
            s->strm->total_in_lo32++;
            if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
         }
         ;
         for(i = 0; i < alphaSize; i++) {
            while(((Bool) 1)) {
               if(curr < 1 || curr > 20) {
                  retVal = (-4);
                  goto save_state_and_return;
               }
               ;
               ;
               case 34:
               s->state = 34;
               while(((Bool) 1)) {
                  if(s->bsLive >= 1) {
                     UInt32 v;
                     v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                     s->bsLive -= 1;
                     uc = v;
                     break;
                  }
                  if(s->strm->avail_in == 0) {
                     retVal = 0;
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                  s->bsLive += 8;
                  s->strm->next_in++;
                  s->strm->avail_in--;
                  s->strm->total_in_lo32++;
                  if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
               }
               ;
               if(uc == 0) break;
               case 35:
               s->state = 35;
               while(((Bool) 1)) {
                  if(s->bsLive >= 1) {
                     UInt32 v;
                     v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                     s->bsLive -= 1;
                     uc = v;
                     break;
                  }
                  if(s->strm->avail_in == 0) {
                     retVal = 0;
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                  s->bsLive += 8;
                  s->strm->next_in++;
                  s->strm->avail_in--;
                  s->strm->total_in_lo32++;
                  if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
               }
               ;
               if(uc == 0) curr++;
               else curr--;
            }
            s->len[t][i] = curr;
         }
      }
      /*--- Create the Huffman decoding tables ---*/
      for(t = 0; t < nGroups; t++) {
         minLen = 32;
         maxLen = 0;
         for(i = 0; i < alphaSize; i++) {
            if(s->len[t][i] > maxLen) maxLen = s->len[t][i];
            if(s->len[t][i] < minLen) minLen = s->len[t][i];
         }
         clava_dcg_global[ 37 ]++;
         BZ2_hbCreateDecodeTables(&(s->limit[t][0]), &(s->base[t][0]), &(s->perm[t][0]), &(s->len[t][0]), minLen, maxLen, alphaSize);
         s->minLens[t] = minLen;
      }
      /*--- Now the MTF values ---*/
      EOB = s->nInUse + 1;
      nblockMAX = 100000 * s->blockSize100k;
      groupNo = -1;
      groupPos = 0;
      for(i = 0; i <= 255; i++) s->unzftab[i] = 0;
      /*-- MTF init --*/
      {
         Int32 ii, jj, kk;
         kk = 4096 - 1;
         for(ii = 256 / 16 - 1; ii >= 0; ii--) {
            for(jj = 16 - 1; jj >= 0; jj--) {
               s->mtfa[kk] = (UChar) (ii * 16 + jj);
               kk--;
            }
            s->mtfbase[ii] = kk + 1;
         }
      }
      /*-- end MTF init --*/
      nblock = 0;
      {
         if(groupPos == 0) {
            groupNo++;
            if(groupNo >= nSelectors) {
               retVal = (-4);
               goto save_state_and_return;
            }
            ;
            ;
            groupPos = 50;
            gSel = s->selector[groupNo];
            gMinlen = s->minLens[gSel];
            gLimit = &(s->limit[gSel][0]);
            gPerm = &(s->perm[gSel][0]);
            gBase = &(s->base[gSel][0]);
         }
         groupPos--;
         zn = gMinlen;
         case 36:
         s->state = 36;
         while(((Bool) 1)) {
            if(s->bsLive >= zn) {
               UInt32 v;
               v = (s->bsBuff >> (s->bsLive - zn)) & ((1 << zn) - 1);
               s->bsLive -= zn;
               zvec = v;
               break;
            }
            if(s->strm->avail_in == 0) {
               retVal = 0;
               goto save_state_and_return;
            }
            ;
            ;
            s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
            s->bsLive += 8;
            s->strm->next_in++;
            s->strm->avail_in--;
            s->strm->total_in_lo32++;
            if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
         }
         ;
         while(1) {
            if(zn > 20) {
               retVal = (-4);
               goto save_state_and_return;
            }
            ;
            ;
            if(zvec <= gLimit[zn]) break;
            zn++;
            case 37:
            s->state = 37;
            while(((Bool) 1)) {
               if(s->bsLive >= 1) {
                  UInt32 v;
                  v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                  s->bsLive -= 1;
                  zj = v;
                  break;
               }
               if(s->strm->avail_in == 0) {
                  retVal = 0;
                  goto save_state_and_return;
               }
               ;
               ;
               s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
               s->bsLive += 8;
               s->strm->next_in++;
               s->strm->avail_in--;
               s->strm->total_in_lo32++;
               if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
            }
            ;
            zvec = (zvec << 1) | zj;
         }
         ;
         if(zvec - gBase[zn] < 0 || zvec - gBase[zn] >= 258) {
            retVal = (-4);
            goto save_state_and_return;
         }
         ;
         ;
         nextSym = gPerm[zvec - gBase[zn]];
      }
      ;
      while(((Bool) 1)) {
         if(nextSym == EOB) break;
         if(nextSym == 0 || nextSym == 1) {
            es = -1;
            N = 1;
            do  {
               if(nextSym == 0) es = es + (0 + 1) * N;
               else if(nextSym == 1) es = es + (1 + 1) * N;
               N = N * 2;
               {
                  if(groupPos == 0) {
                     groupNo++;
                     if(groupNo >= nSelectors) {
                        retVal = (-4);
                        goto save_state_and_return;
                     }
                     ;
                     ;
                     groupPos = 50;
                     gSel = s->selector[groupNo];
                     gMinlen = s->minLens[gSel];
                     gLimit = &(s->limit[gSel][0]);
                     gPerm = &(s->perm[gSel][0]);
                     gBase = &(s->base[gSel][0]);
                  }
                  groupPos--;
                  zn = gMinlen;
                  case 38:
                  s->state = 38;
                  while(((Bool) 1)) {
                     if(s->bsLive >= zn) {
                        UInt32 v;
                        v = (s->bsBuff >> (s->bsLive - zn)) & ((1 << zn) - 1);
                        s->bsLive -= zn;
                        zvec = v;
                        break;
                     }
                     if(s->strm->avail_in == 0) {
                        retVal = 0;
                        goto save_state_and_return;
                     }
                     ;
                     ;
                     s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                     s->bsLive += 8;
                     s->strm->next_in++;
                     s->strm->avail_in--;
                     s->strm->total_in_lo32++;
                     if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
                  }
                  ;
                  while(1) {
                     if(zn > 20) {
                        retVal = (-4);
                        goto save_state_and_return;
                     }
                     ;
                     ;
                     if(zvec <= gLimit[zn]) break;
                     zn++;
                     case 39:
                     s->state = 39;
                     while(((Bool) 1)) {
                        if(s->bsLive >= 1) {
                           UInt32 v;
                           v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                           s->bsLive -= 1;
                           zj = v;
                           break;
                        }
                        if(s->strm->avail_in == 0) {
                           retVal = 0;
                           goto save_state_and_return;
                        }
                        ;
                        ;
                        s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                        s->bsLive += 8;
                        s->strm->next_in++;
                        s->strm->avail_in--;
                        s->strm->total_in_lo32++;
                        if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
                     }
                     ;
                     zvec = (zvec << 1) | zj;
                  }
                  ;
                  if(zvec - gBase[zn] < 0 || zvec - gBase[zn] >= 258) {
                     retVal = (-4);
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  nextSym = gPerm[zvec - gBase[zn]];
               }
               ;
            }
            while (nextSym == 0 || nextSym == 1);
            es++;
            uc = s->seqToUnseq[s->mtfa[s->mtfbase[0]]];
            s->unzftab[uc] += es;
            if(s->smallDecompress) while(es > 0) {
               if(nblock >= nblockMAX) {
                  retVal = (-4);
                  goto save_state_and_return;
               }
               ;
               ;
               s->ll16[nblock] = (UInt16) uc;
               nblock++;
               es--;
            }
            else while(es > 0) {
               if(nblock >= nblockMAX) {
                  retVal = (-4);
                  goto save_state_and_return;
               }
               ;
               ;
               s->tt[nblock] = (UInt32) uc;
               nblock++;
               es--;
            }
            ;
            continue;
         }
         else {
            if(nblock >= nblockMAX) {
               retVal = (-4);
               goto save_state_and_return;
            }
            ;
            ;
            /*-- uc = MTF ( nextSym-1 ) --*/
            {
               Int32 ii, jj, kk, pp, lno, off;
               UInt32 nn;
               nn = (UInt32) (nextSym - 1);
               if(nn < 16) {
                  /*avoid general-case expense*/
                  pp = s->mtfbase[0];
                  uc = s->mtfa[pp + nn];
                  while(nn > 3) {
                     Int32 z = pp + nn;
                     s->mtfa[(z)] = s->mtfa[(z) - 1];
                     s->mtfa[(z) - 1] = s->mtfa[(z) - 2];
                     s->mtfa[(z) - 2] = s->mtfa[(z) - 3];
                     s->mtfa[(z) - 3] = s->mtfa[(z) - 4];
                     nn -= 4;
                  }
                  while(nn > 0) {
                     s->mtfa[(pp + nn)] = s->mtfa[(pp + nn) - 1];
                     nn--;
                  }
                  ;
                  s->mtfa[pp] = uc;
               }
               else {
                  /*general case*/
                  lno = nn / 16;
                  off = nn % 16;
                  pp = s->mtfbase[lno] + off;
                  uc = s->mtfa[pp];
                  while(pp > s->mtfbase[lno]) {
                     s->mtfa[pp] = s->mtfa[pp - 1];
                     pp--;
                  }
                  ;
                  s->mtfbase[lno]++;
                  while(lno > 0) {
                     s->mtfbase[lno]--;
                     s->mtfa[s->mtfbase[lno]] = s->mtfa[s->mtfbase[lno - 1] + 16 - 1];
                     lno--;
                  }
                  s->mtfbase[0]--;
                  s->mtfa[s->mtfbase[0]] = uc;
                  if(s->mtfbase[0] == 0) {
                     kk = 4096 - 1;
                     for(ii = 256 / 16 - 1; ii >= 0; ii--) {
                        for(jj = 16 - 1; jj >= 0; jj--) {
                           s->mtfa[kk] = s->mtfa[s->mtfbase[ii] + jj];
                           kk--;
                        }
                        s->mtfbase[ii] = kk + 1;
                     }
                  }
               }
            }
            /*-- end uc = MTF ( nextSym-1 ) --*/
            s->unzftab[s->seqToUnseq[uc]]++;
            if(s->smallDecompress) s->ll16[nblock] = (UInt16) (s->seqToUnseq[uc]);
            else s->tt[nblock] = (UInt32) (s->seqToUnseq[uc]);
            nblock++;
            {
               if(groupPos == 0) {
                  groupNo++;
                  if(groupNo >= nSelectors) {
                     retVal = (-4);
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  groupPos = 50;
                  gSel = s->selector[groupNo];
                  gMinlen = s->minLens[gSel];
                  gLimit = &(s->limit[gSel][0]);
                  gPerm = &(s->perm[gSel][0]);
                  gBase = &(s->base[gSel][0]);
               }
               groupPos--;
               zn = gMinlen;
               case 40:
               s->state = 40;
               while(((Bool) 1)) {
                  if(s->bsLive >= zn) {
                     UInt32 v;
                     v = (s->bsBuff >> (s->bsLive - zn)) & ((1 << zn) - 1);
                     s->bsLive -= zn;
                     zvec = v;
                     break;
                  }
                  if(s->strm->avail_in == 0) {
                     retVal = 0;
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                  s->bsLive += 8;
                  s->strm->next_in++;
                  s->strm->avail_in--;
                  s->strm->total_in_lo32++;
                  if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
               }
               ;
               while(1) {
                  if(zn > 20) {
                     retVal = (-4);
                     goto save_state_and_return;
                  }
                  ;
                  ;
                  if(zvec <= gLimit[zn]) break;
                  zn++;
                  case 41:
                  s->state = 41;
                  while(((Bool) 1)) {
                     if(s->bsLive >= 1) {
                        UInt32 v;
                        v = (s->bsBuff >> (s->bsLive - 1)) & ((1 << 1) - 1);
                        s->bsLive -= 1;
                        zj = v;
                        break;
                     }
                     if(s->strm->avail_in == 0) {
                        retVal = 0;
                        goto save_state_and_return;
                     }
                     ;
                     ;
                     s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
                     s->bsLive += 8;
                     s->strm->next_in++;
                     s->strm->avail_in--;
                     s->strm->total_in_lo32++;
                     if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
                  }
                  ;
                  zvec = (zvec << 1) | zj;
               }
               ;
               if(zvec - gBase[zn] < 0 || zvec - gBase[zn] >= 258) {
                  retVal = (-4);
                  goto save_state_and_return;
               }
               ;
               ;
               nextSym = gPerm[zvec - gBase[zn]];
            }
            ;
            continue;
         }
      }
      /*Now we know what nblock is, we can do a better sanity
      check on s->origPtr.
      */
      if(s->origPtr < 0 || s->origPtr >= nblock) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      s->state_out_len = 0;
      s->state_out_ch = 0;
      {
         s->calculatedBlockCRC = 0xffffffffL;
      }
      ;
      s->state = 2;
      if(s->verbosity >= 2) {
         clava_dcg_global[ 35 ]++;
         fprintf(stderr, "rt+rld");
      }
      /*-- Set up cftab to facilitate generation of T^(-1) --*/
      s->cftab[0] = 0;
      for(i = 1; i <= 256; i++) s->cftab[i] = s->unzftab[i - 1];
      for(i = 1; i <= 256; i++) s->cftab[i] += s->cftab[i - 1];
      if(s->smallDecompress) {
         /*-- Make a copy of cftab, used in generation of T --*/
         for(i = 0; i <= 256; i++) s->cftabCopy[i] = s->cftab[i];
         /*-- compute the T vector --*/
         for(i = 0; i < nblock; i++) {
            uc = (UChar) (s->ll16[i]);
            {
               s->ll16[i] = (UInt16) (s->cftabCopy[uc] & 0x0000ffff);
               {
                  if(((i) & 0x1) == 0) s->ll4[(i) >> 1] = (s->ll4[(i) >> 1] & 0xf0) | (s->cftabCopy[uc] >> 16);
                  else s->ll4[(i) >> 1] = (s->ll4[(i) >> 1] & 0x0f) | ((s->cftabCopy[uc] >> 16) << 4);
               }
               ;
            }
            ;
            s->cftabCopy[uc]++;
         }
         /*-- Compute T^(-1) by pointer reversal on T --*/
         i = s->origPtr;
         j = (((UInt32) s->ll16[i]) | (((((UInt32) (s->ll4[(i) >> 1])) >> (((i) << 2) & 0x4)) & 0xF) << 16));
         do  {
            Int32 tmp = (((UInt32) s->ll16[j]) | (((((UInt32) (s->ll4[(j) >> 1])) >> (((j) << 2) & 0x4)) & 0xF) << 16));
            {
               s->ll16[j] = (UInt16) (i & 0x0000ffff);
               {
                  if(((j) & 0x1) == 0) s->ll4[(j) >> 1] = (s->ll4[(j) >> 1] & 0xf0) | (i >> 16);
                  else s->ll4[(j) >> 1] = (s->ll4[(j) >> 1] & 0x0f) | ((i >> 16) << 4);
               }
               ;
            }
            ;
            i = j;
            j = tmp;
         }
         while (i != s->origPtr);
         s->tPos = s->origPtr;
         s->nblock_used = 0;
         if(s->blockRandomised) {
            s->rNToGo = 0;
            s->rTPos = 0;
            clava_dcg_global[ 38 ]++;
            s->k0 = BZ2_indexIntoF(s->tPos, s->cftab);
            s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
            ;
            s->nblock_used++;
            if(s->rNToGo == 0) {
               s->rNToGo = BZ2_rNums[s->rTPos];
               s->rTPos++;
               if(s->rTPos == 512) s->rTPos = 0;
            }
            s->rNToGo--;
            ;
            s->k0 ^= ((s->rNToGo == 1) ? 1 : 0);
         }
         else {
            clava_dcg_global[ 38 ]++;
            s->k0 = BZ2_indexIntoF(s->tPos, s->cftab);
            s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
            ;
            s->nblock_used++;
         }
      }
      else {
         /*-- compute the T^(-1) vector --*/
         for(i = 0; i < nblock; i++) {
            uc = (UChar) (s->tt[i] & 0xff);
            s->tt[s->cftab[uc]] |= (i << 8);
            s->cftab[uc]++;
         }
         s->tPos = s->tt[s->origPtr] >> 8;
         s->nblock_used = 0;
         if(s->blockRandomised) {
            s->rNToGo = 0;
            s->rTPos = 0;
            s->tPos = s->tt[s->tPos];
            s->k0 = (UChar) (s->tPos & 0xff);
            s->tPos >>= 8;
            ;
            s->nblock_used++;
            if(s->rNToGo == 0) {
               s->rNToGo = BZ2_rNums[s->rTPos];
               s->rTPos++;
               if(s->rTPos == 512) s->rTPos = 0;
            }
            s->rNToGo--;
            ;
            s->k0 ^= ((s->rNToGo == 1) ? 1 : 0);
         }
         else {
            s->tPos = s->tt[s->tPos];
            s->k0 = (UChar) (s->tPos & 0xff);
            s->tPos >>= 8;
            ;
            s->nblock_used++;
         }
      }
      {
         retVal = 0;
         goto save_state_and_return;
      }
      ;
      ;
      endhdr_2:
      case 42:
      s->state = 42;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x72) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 43:
      s->state = 43;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x45) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 44:
      s->state = 44;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x38) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 45:
      s->state = 45;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x50) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      case 46:
      s->state = 46;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      if(uc != 0x90) {
         retVal = (-4);
         goto save_state_and_return;
      }
      ;
      ;
      s->storedCombinedCRC = 0;
      case 47:
      s->state = 47;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32) uc);
      case 48:
      s->state = 48;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32) uc);
      case 49:
      s->state = 49;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32) uc);
      case 50:
      s->state = 50;
      while(((Bool) 1)) {
         if(s->bsLive >= 8) {
            UInt32 v;
            v = (s->bsBuff >> (s->bsLive - 8)) & ((1 << 8) - 1);
            s->bsLive -= 8;
            uc = v;
            break;
         }
         if(s->strm->avail_in == 0) {
            retVal = 0;
            goto save_state_and_return;
         }
         ;
         ;
         s->bsBuff = (s->bsBuff << 8) | ((UInt32) (*((UChar *) (s->strm->next_in))));
         s->bsLive += 8;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
      ;
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32) uc);
      s->state = 1;
      {
         retVal = 4;
         goto save_state_and_return;
      }
      ;
      ;
      default:
      {
         if(!(((Bool) 0))) {
            clava_dcg_global[ 39 ]++;
            BZ2_bz__AssertH__fail(4001);
         }
      }
      ;
   }
   {
      if(!(((Bool) 0))) {
         clava_dcg_global[ 39 ]++;
         BZ2_bz__AssertH__fail(4002);
      }
   }
   ;
   save_state_and_return:
   s->save_i = i;
   s->save_j = j;
   s->save_t = t;
   s->save_alphaSize = alphaSize;
   s->save_nGroups = nGroups;
   s->save_nSelectors = nSelectors;
   s->save_EOB = EOB;
   s->save_groupNo = groupNo;
   s->save_groupPos = groupPos;
   s->save_nextSym = nextSym;
   s->save_nblockMAX = nblockMAX;
   s->save_nblock = nblock;
   s->save_es = es;
   s->save_N = N;
   s->save_curr = curr;
   s->save_zt = zt;
   s->save_zn = zn;
   s->save_zvec = zvec;
   s->save_zj = zj;
   s->save_gSel = gSel;
   s->save_gMinlen = gMinlen;
   s->save_gLimit = gLimit;
   s->save_gBase = gBase;
   s->save_gPerm = gPerm;
   
   return retVal;
}

/*-------------------------------------------------------------*/
/*--- end                                      decompress.c ---*/
/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/
/*--- Library top-level functions.                          ---*/
/*---                                               bzlib.c ---*/
/*-------------------------------------------------------------*/
/*---------------------------------------------------*/
/*--- Compression stuff                           ---*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
void BZ2_bz__AssertH__fail(int errcode) {
   clava_dcg_global[ 40 ]++;
   clava_dcg_global[ 41 ]++;
   fprintf(stderr, "\n\nbzip2/libbzip2: internal error number %d.\nThis is a bug in bzip2/libbzip2, %s.\nPlease report it to me at: jseward@acm.org.  If this happened\nwhen you were using some program which uses libbzip2 as a\ncomponent, you should also report this bug to the author(s)\nof that program.  Please make an effort to report this bug;\ntimely and accurate bug reports eventually lead to higher\nquality software.  Thanks.  Julian Seward, 30 December 2001.\n\n", errcode, BZ2_bzlibVersion());
   if(errcode == 1007) {
      clava_dcg_global[ 40 ]++;
      fprintf(stderr, "\n*** A special note about internal error number 1007 ***\n\nExperience suggests that a common cause of i.e. 1007\nis unreliable memory or other hardware.  The 1007 assertion\njust happens to cross-check the results of huge numbers of\nmemory reads/writes, and so acts (unintendedly) as a stress\ntest of your memory system.\n\nI suggest the following: try compressing the file again,\npossibly monitoring progress in detail with the -vv flag.\n\n* If the error cannot be reproduced, and/or happens at different\n  points in compression, you may have a flaky memory system.\n  Try a memory-test program.  I have used Memtest86\n  (www.memtest86.com).  At the time of writing it is free (GPLd).\n  Memtest86 tests memory much more thorougly than your BIOSs\n  power-on test, and may find failures that the BIOS doesn't.\n\n* If the error can be repeatably reproduced, this is a bug in\n  bzip2, and I would very much like to hear about it.  Please\n  let me know, and, ideally, save a copy of the file causing the\n  problem -- without which I will be unable to investigate it.\n\n");
   }
   clava_dcg_global[ 42 ]++;
   exit(3);
}

/*---------------------------------------------------*/
static int bz_config_ok() {
   if(sizeof(int) != 4) 
   return 0;
   if(sizeof(short) != 2) 
   return 0;
   if(sizeof(char) != 1) 
   return 0;
   
   return 1;
}

/*---------------------------------------------------*/
static void * default_bzalloc(void *opaque, Int32 items, Int32 size) {
   clava_dcg_global[ 43 ]++;
   void *v = malloc(items * size);
   
   return v;
}

static void default_bzfree(void *opaque, void *addr) {
   if(addr != ((void *) 0)) {
      clava_dcg_global[ 44 ]++;
      free(addr);
   }
}

/*---------------------------------------------------*/
static void prepare_new_block(EState *s) {
   Int32 i;
   s->nblock = 0;
   s->numZ = 0;
   s->state_out_pos = 0;
   {
      s->blockCRC = 0xffffffffL;
   }
   ;
   for(i = 0; i < 256; i++) s->inUse[i] = ((Bool) 0);
   s->blockNo++;
}

/*---------------------------------------------------*/
static void init_RL(EState *s) {
   s->state_in_ch = 256;
   s->state_in_len = 0;
}

static Bool isempty_RL(EState *s) {
   if(s->state_in_ch < 256 && s->state_in_len > 0) 
   return ((Bool) 0);
   else 
   return ((Bool) 1);
}

/*---------------------------------------------------*/
int BZ2_bzCompressInit(bz_stream *strm, int blockSize100k, int verbosity, int workFactor) {
   Int32 n;
   EState *s;
   clava_dcg_global[ 45 ]++;
   if(!bz_config_ok()) 
   return (-9);
   if(strm == ((void *) 0) || blockSize100k < 1 || blockSize100k > 9 || workFactor < 0 || workFactor > 250) 
   return (-2);
   if(workFactor == 0) workFactor = 30;
   if(strm->bzalloc == ((void *) 0)) strm->bzalloc = default_bzalloc;
   if(strm->bzfree == ((void *) 0)) strm->bzfree = default_bzfree;
   clava_dcg_global[ 46 ]++;
   s = (strm->bzalloc)(strm->opaque, (sizeof(EState)), 1);
   if(s == ((void *) 0)) 
   return (-3);
   s->strm = strm;
   s->arr1 = ((void *) 0);
   s->arr2 = ((void *) 0);
   s->ftab = ((void *) 0);
   n = 100000 * blockSize100k;
   clava_dcg_global[ 46 ]++;
   s->arr1 = (strm->bzalloc)(strm->opaque, (n * sizeof(UInt32)), 1);
   clava_dcg_global[ 46 ]++;
   s->arr2 = (strm->bzalloc)(strm->opaque, ((n + (2 + 12 + 18 + 2)) * sizeof(UInt32)), 1);
   clava_dcg_global[ 46 ]++;
   s->ftab = (strm->bzalloc)(strm->opaque, (65537 * sizeof(UInt32)), 1);
   if(s->arr1 == ((void *) 0) || s->arr2 == ((void *) 0) || s->ftab == ((void *) 0)) {
      if(s->arr1 != ((void *) 0)) {
         clava_dcg_global[ 46 ]++;
         (strm->bzfree)(strm->opaque, (s->arr1));
      }
      if(s->arr2 != ((void *) 0)) {
         clava_dcg_global[ 46 ]++;
         (strm->bzfree)(strm->opaque, (s->arr2));
      }
      if(s->ftab != ((void *) 0)) {
         clava_dcg_global[ 46 ]++;
         (strm->bzfree)(strm->opaque, (s->ftab));
      }
      if(s != ((void *) 0)) {
         clava_dcg_global[ 46 ]++;
         (strm->bzfree)(strm->opaque, (s));
      }
      
      return (-3);
   }
   s->blockNo = 0;
   s->state = 2;
   s->mode = 2;
   s->combinedCRC = 0;
   s->blockSize100k = blockSize100k;
   s->nblockMAX = 100000 * blockSize100k - 19;
   s->verbosity = verbosity;
   s->workFactor = workFactor;
   s->block = (UChar *) s->arr2;
   s->mtfv = (UInt16 *) s->arr1;
   s->zbits = ((void *) 0);
   s->ptr = (UInt32 *) s->arr1;
   strm->state = s;
   strm->total_in_lo32 = 0;
   strm->total_in_hi32 = 0;
   strm->total_out_lo32 = 0;
   strm->total_out_hi32 = 0;
   clava_dcg_global[ 47 ]++;
   init_RL(s);
   clava_dcg_global[ 48 ]++;
   prepare_new_block(s);
   
   return 0;
}

/*---------------------------------------------------*/
static void add_pair_to_block(EState *s) {
   Int32 i;
   UChar ch = (UChar) (s->state_in_ch);
   for(i = 0; i < s->state_in_len; i++) {
      {
         s->blockCRC = (s->blockCRC << 8) ^ BZ2_crc32Table[(s->blockCRC >> 24) ^ ((UChar) ch)];
      }
      ;
   }
   s->inUse[s->state_in_ch] = ((Bool) 1);
   switch (s->state_in_len) {
      case 1:
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      break;
      case 2:
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      break;
      case 3:
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      break;
      default:
      s->inUse[s->state_in_len - 4] = ((Bool) 1);
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = (UChar) ch;
      s->nblock++;
      s->block[s->nblock] = ((UChar) (s->state_in_len - 4));
      s->nblock++;
      break;
   }
}

/*---------------------------------------------------*/
static void flush_RL(EState *s) {
   if(s->state_in_ch < 256) {
      clava_dcg_global[ 49 ]++;
      add_pair_to_block(s);
   }
   clava_dcg_global[ 50 ]++;
   init_RL(s);
}

/*---------------------------------------------------*/
/*-- fast track the common case --*/
/*-- general, uncommon cases --*/
/*---------------------------------------------------*/
static Bool copy_input_until_stop(EState *s) {
   Bool progress_in = ((Bool) 0);
   if(s->mode == 2) {
      /*-- fast track the common case --*/
      while(((Bool) 1)) {
         /*-- block full? --*/
         if(s->nblock >= s->nblockMAX) break;
         /*-- no input? --*/
         if(s->strm->avail_in == 0) break;
         progress_in = ((Bool) 1);
         {
            UInt32 zchh = (UInt32) ((UInt32) (*((UChar *) (s->strm->next_in))));
            if(zchh != s->state_in_ch && s->state_in_len == 1) {
               UChar ch = (UChar) (s->state_in_ch);
               {
                  s->blockCRC = (s->blockCRC << 8) ^ BZ2_crc32Table[(s->blockCRC >> 24) ^ ((UChar) ch)];
               }
               ;
               s->inUse[s->state_in_ch] = ((Bool) 1);
               s->block[s->nblock] = (UChar) ch;
               s->nblock++;
               s->state_in_ch = zchh;
            }
            else if(zchh != s->state_in_ch || s->state_in_len == 255) {
               if(s->state_in_ch < 256) {
                  clava_dcg_global[ 51 ]++;
                  add_pair_to_block(s);
               }
               s->state_in_ch = zchh;
               s->state_in_len = 1;
            }
            else {
               s->state_in_len++;
            }
         }
         ;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
      }
   }
   else {
      /*-- general, uncommon case --*/
      while(((Bool) 1)) {
         /*-- block full? --*/
         if(s->nblock >= s->nblockMAX) break;
         /*-- no input? --*/
         if(s->strm->avail_in == 0) break;
         /*-- flush/finish end? --*/
         if(s->avail_in_expect == 0) break;
         progress_in = ((Bool) 1);
         {
            UInt32 zchh = (UInt32) ((UInt32) (*((UChar *) (s->strm->next_in))));
            if(zchh != s->state_in_ch && s->state_in_len == 1) {
               UChar ch = (UChar) (s->state_in_ch);
               {
                  s->blockCRC = (s->blockCRC << 8) ^ BZ2_crc32Table[(s->blockCRC >> 24) ^ ((UChar) ch)];
               }
               ;
               s->inUse[s->state_in_ch] = ((Bool) 1);
               s->block[s->nblock] = (UChar) ch;
               s->nblock++;
               s->state_in_ch = zchh;
            }
            else if(zchh != s->state_in_ch || s->state_in_len == 255) {
               if(s->state_in_ch < 256) {
                  clava_dcg_global[ 51 ]++;
                  add_pair_to_block(s);
               }
               s->state_in_ch = zchh;
               s->state_in_len = 1;
            }
            else {
               s->state_in_len++;
            }
         }
         ;
         s->strm->next_in++;
         s->strm->avail_in--;
         s->strm->total_in_lo32++;
         if(s->strm->total_in_lo32 == 0) s->strm->total_in_hi32++;
         s->avail_in_expect--;
      }
   }
   
   return progress_in;
}

/*---------------------------------------------------*/
static Bool copy_output_until_stop(EState *s) {
   Bool progress_out = ((Bool) 0);
   while(((Bool) 1)) {
      /*-- no output space? --*/
      if(s->strm->avail_out == 0) break;
      /*-- block done? --*/
      if(s->state_out_pos >= s->numZ) break;
      progress_out = ((Bool) 1);
      *(s->strm->next_out) = s->zbits[s->state_out_pos];
      s->state_out_pos++;
      s->strm->avail_out--;
      s->strm->next_out++;
      s->strm->total_out_lo32++;
      if(s->strm->total_out_lo32 == 0) s->strm->total_out_hi32++;
   }
   
   return progress_out;
}

/*---------------------------------------------------*/
static Bool handle_compress(bz_stream *strm) {
   Bool progress_in = ((Bool) 0);
   Bool progress_out = ((Bool) 0);
   EState *s = strm->state;
   while(((Bool) 1)) {
      if(s->state == 1) {
         clava_dcg_global[ 52 ]++;
         progress_out |= copy_output_until_stop(s);
         if(s->state_out_pos < s->numZ) break;
         clava_dcg_global[ 53 ]++;
         if(s->mode == 4 && s->avail_in_expect == 0 && isempty_RL(s)) break;
         clava_dcg_global[ 54 ]++;
         prepare_new_block(s);
         s->state = 2;
         clava_dcg_global[ 53 ]++;
         if(s->mode == 3 && s->avail_in_expect == 0 && isempty_RL(s)) break;
      }
      if(s->state == 2) {
         clava_dcg_global[ 55 ]++;
         progress_in |= copy_input_until_stop(s);
         if(s->mode != 2 && s->avail_in_expect == 0) {
            clava_dcg_global[ 56 ]++;
            flush_RL(s);
            clava_dcg_global[ 57 ]++;
            BZ2_compressBlock(s, (Bool) (s->mode == 4));
            s->state = 1;
         }
         else if(s->nblock >= s->nblockMAX) {
            clava_dcg_global[ 57 ]++;
            BZ2_compressBlock(s, ((Bool) 0));
            s->state = 1;
         }
         else if(s->strm->avail_in == 0) {
            break;
         }
      }
   }
   
   return progress_in || progress_out;
}

/*---------------------------------------------------*/
int BZ2_bzCompress(bz_stream *strm, int action) {
   Bool progress;
   EState *s;
   if(strm == ((void *) 0)) 
   return (-2);
   s = strm->state;
   if(s == ((void *) 0)) 
   return (-2);
   if(s->strm != strm) 
   return (-2);
   preswitch:
   switch (s->mode) {
      case 1:
      
      return (-1);
      case 2:
      if(action == 0) {
         clava_dcg_global[ 58 ]++;
         progress = handle_compress(strm);
         
         return progress ? 1 : (-2);
      }
      else if(action == 1) {
         s->avail_in_expect = strm->avail_in;
         s->mode = 3;
         goto preswitch;
      }
      else if(action == 2) {
         s->avail_in_expect = strm->avail_in;
         s->mode = 4;
         goto preswitch;
      }
      else 
      return (-2);
      case 3:
      if(action != 1) 
      return (-1);
      if(s->avail_in_expect != s->strm->avail_in) 
      return (-1);
      clava_dcg_global[ 58 ]++;
      progress = handle_compress(strm);
      clava_dcg_global[ 59 ]++;
      if(s->avail_in_expect > 0 || !isempty_RL(s) || s->state_out_pos < s->numZ) 
      return 2;
      s->mode = 2;
      
      return 1;
      case 4:
      if(action != 2) 
      return (-1);
      if(s->avail_in_expect != s->strm->avail_in) 
      return (-1);
      clava_dcg_global[ 58 ]++;
      progress = handle_compress(strm);
      if(!progress) 
      return (-1);
      clava_dcg_global[ 59 ]++;
      if(s->avail_in_expect > 0 || !isempty_RL(s) || s->state_out_pos < s->numZ) 
      return 3;
      s->mode = 1;
      
      return 4;
   }
   
   return 0;
   /*--not reached--*/
}

/*---------------------------------------------------*/
int BZ2_bzCompressEnd(bz_stream *strm) {
   EState *s;
   if(strm == ((void *) 0)) 
   return (-2);
   s = strm->state;
   if(s == ((void *) 0)) 
   return (-2);
   if(s->strm != strm) 
   return (-2);
   if(s->arr1 != ((void *) 0)) {
      clava_dcg_global[ 60 ]++;
      (strm->bzfree)(strm->opaque, (s->arr1));
   }
   if(s->arr2 != ((void *) 0)) {
      clava_dcg_global[ 60 ]++;
      (strm->bzfree)(strm->opaque, (s->arr2));
   }
   if(s->ftab != ((void *) 0)) {
      clava_dcg_global[ 60 ]++;
      (strm->bzfree)(strm->opaque, (s->ftab));
   }
   clava_dcg_global[ 60 ]++;
   (strm->bzfree)(strm->opaque, (strm->state));
   strm->state = ((void *) 0);
   
   return 0;
}

/*---------------------------------------------------*/
/*--- Decompression stuff                         ---*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
int BZ2_bzDecompressInit(bz_stream *strm, int verbosity, int small) {
   DState *s;
   clava_dcg_global[ 61 ]++;
   if(!bz_config_ok()) 
   return (-9);
   if(strm == ((void *) 0)) 
   return (-2);
   if(small != 0 && small != 1) 
   return (-2);
   if(verbosity < 0 || verbosity > 4) 
   return (-2);
   if(strm->bzalloc == ((void *) 0)) strm->bzalloc = default_bzalloc;
   if(strm->bzfree == ((void *) 0)) strm->bzfree = default_bzfree;
   clava_dcg_global[ 62 ]++;
   s = (strm->bzalloc)(strm->opaque, (sizeof(DState)), 1);
   if(s == ((void *) 0)) 
   return (-3);
   s->strm = strm;
   strm->state = s;
   s->state = 10;
   s->bsLive = 0;
   s->bsBuff = 0;
   s->calculatedCombinedCRC = 0;
   strm->total_in_lo32 = 0;
   strm->total_in_hi32 = 0;
   strm->total_out_lo32 = 0;
   strm->total_out_hi32 = 0;
   s->smallDecompress = (Bool) small;
   s->ll4 = ((void *) 0);
   s->ll16 = ((void *) 0);
   s->tt = ((void *) 0);
   s->currBlockNo = 0;
   s->verbosity = verbosity;
   
   return 0;
}

/*---------------------------------------------------*/
static void unRLE_obuf_to_output_FAST(DState *s) {
   UChar k1;
   if(s->blockRandomised) {
      while(((Bool) 1)) {
         /*try to finish existing run*/
         while(((Bool) 1)) {
            if(s->strm->avail_out == 0) 
            return;
            if(s->state_out_len == 0) break;
            *((UChar *) (s->strm->next_out)) = s->state_out_ch;
            {
               s->calculatedBlockCRC = (s->calculatedBlockCRC << 8) ^ BZ2_crc32Table[(s->calculatedBlockCRC >> 24) ^ ((UChar) s->state_out_ch)];
            }
            ;
            s->state_out_len--;
            s->strm->next_out++;
            s->strm->avail_out--;
            s->strm->total_out_lo32++;
            if(s->strm->total_out_lo32 == 0) s->strm->total_out_hi32++;
         }
         /*can a new run be started?*/
         if(s->nblock_used == s->save_nblock + 1) 
         return;
         s->state_out_len = 1;
         s->state_out_ch = s->k0;
         s->tPos = s->tt[s->tPos];
         k1 = (UChar) (s->tPos & 0xff);
         s->tPos >>= 8;
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 2;
         s->tPos = s->tt[s->tPos];
         k1 = (UChar) (s->tPos & 0xff);
         s->tPos >>= 8;
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 3;
         s->tPos = s->tt[s->tPos];
         k1 = (UChar) (s->tPos & 0xff);
         s->tPos >>= 8;
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->tPos = s->tt[s->tPos];
         k1 = (UChar) (s->tPos & 0xff);
         s->tPos >>= 8;
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         s->state_out_len = ((Int32) k1) + 4;
         s->tPos = s->tt[s->tPos];
         s->k0 = (UChar) (s->tPos & 0xff);
         s->tPos >>= 8;
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         s->k0 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
      }
   }
   else {
      /*restore*/
      UInt32 c_calculatedBlockCRC = s->calculatedBlockCRC;
      UChar c_state_out_ch = s->state_out_ch;
      Int32 c_state_out_len = s->state_out_len;
      Int32 c_nblock_used = s->nblock_used;
      Int32 c_k0 = s->k0;
      UInt32 *c_tt = s->tt;
      UInt32 c_tPos = s->tPos;
      char *cs_next_out = s->strm->next_out;
      unsigned int cs_avail_out = s->strm->avail_out;
      /*end restore*/
      UInt32 avail_out_INIT = cs_avail_out;
      Int32 s_save_nblockPP = s->save_nblock + 1;
      unsigned int total_out_lo32_old;
      while(((Bool) 1)) {
         /*try to finish existing run*/
         if(c_state_out_len > 0) {
            while(((Bool) 1)) {
               if(cs_avail_out == 0) goto return_notr;
               if(c_state_out_len == 1) break;
               *((UChar *) (cs_next_out)) = c_state_out_ch;
               {
                  c_calculatedBlockCRC = (c_calculatedBlockCRC << 8) ^ BZ2_crc32Table[(c_calculatedBlockCRC >> 24) ^ ((UChar) c_state_out_ch)];
               }
               ;
               c_state_out_len--;
               cs_next_out++;
               cs_avail_out--;
            }
            s_state_out_len_eq_one:
            {
               if(cs_avail_out == 0) {
                  c_state_out_len = 1;
                  goto return_notr;
               }
               ;
               *((UChar *) (cs_next_out)) = c_state_out_ch;
               {
                  c_calculatedBlockCRC = (c_calculatedBlockCRC << 8) ^ BZ2_crc32Table[(c_calculatedBlockCRC >> 24) ^ ((UChar) c_state_out_ch)];
               }
               ;
               cs_next_out++;
               cs_avail_out--;
            }
         }
         /*can a new run be started?*/
         if(c_nblock_used == s_save_nblockPP) {
            c_state_out_len = 0;
            goto return_notr;
         }
         ;
         c_state_out_ch = c_k0;
         c_tPos = c_tt[c_tPos];
         k1 = (UChar) (c_tPos & 0xff);
         c_tPos >>= 8;
         ;
         c_nblock_used++;
         if(k1 != c_k0) {
            c_k0 = k1;
            goto s_state_out_len_eq_one;
         }
         ;
         if(c_nblock_used == s_save_nblockPP) goto s_state_out_len_eq_one;
         c_state_out_len = 2;
         c_tPos = c_tt[c_tPos];
         k1 = (UChar) (c_tPos & 0xff);
         c_tPos >>= 8;
         ;
         c_nblock_used++;
         if(c_nblock_used == s_save_nblockPP) continue;
         if(k1 != c_k0) {
            c_k0 = k1;
            continue;
         }
         ;
         c_state_out_len = 3;
         c_tPos = c_tt[c_tPos];
         k1 = (UChar) (c_tPos & 0xff);
         c_tPos >>= 8;
         ;
         c_nblock_used++;
         if(c_nblock_used == s_save_nblockPP) continue;
         if(k1 != c_k0) {
            c_k0 = k1;
            continue;
         }
         ;
         c_tPos = c_tt[c_tPos];
         k1 = (UChar) (c_tPos & 0xff);
         c_tPos >>= 8;
         ;
         c_nblock_used++;
         c_state_out_len = ((Int32) k1) + 4;
         c_tPos = c_tt[c_tPos];
         c_k0 = (UChar) (c_tPos & 0xff);
         c_tPos >>= 8;
         ;
         c_nblock_used++;
      }
      return_notr:
      total_out_lo32_old = s->strm->total_out_lo32;
      s->strm->total_out_lo32 += (avail_out_INIT - cs_avail_out);
      if(s->strm->total_out_lo32 < total_out_lo32_old) s->strm->total_out_hi32++;
      /*save*/
      s->calculatedBlockCRC = c_calculatedBlockCRC;
      s->state_out_ch = c_state_out_ch;
      s->state_out_len = c_state_out_len;
      s->nblock_used = c_nblock_used;
      s->k0 = c_k0;
      s->tt = c_tt;
      s->tPos = c_tPos;
      s->strm->next_out = cs_next_out;
      s->strm->avail_out = cs_avail_out;
      /*end save*/
   }
}

/*---------------------------------------------------*/
Int32 BZ2_indexIntoF(Int32 indx, Int32 *cftab) {
   Int32 nb, na, mid;
   nb = 0;
   na = 256;
   do  {
      mid = (nb + na) >> 1;
      if(indx >= cftab[mid]) nb = mid;
      else na = mid;
   }
   while (na - nb != 1);
   
   return nb;
}

/*---------------------------------------------------*/
static void unRLE_obuf_to_output_SMALL(DState *s) {
   UChar k1;
   if(s->blockRandomised) {
      while(((Bool) 1)) {
         /*try to finish existing run*/
         while(((Bool) 1)) {
            if(s->strm->avail_out == 0) 
            return;
            if(s->state_out_len == 0) break;
            *((UChar *) (s->strm->next_out)) = s->state_out_ch;
            {
               s->calculatedBlockCRC = (s->calculatedBlockCRC << 8) ^ BZ2_crc32Table[(s->calculatedBlockCRC >> 24) ^ ((UChar) s->state_out_ch)];
            }
            ;
            s->state_out_len--;
            s->strm->next_out++;
            s->strm->avail_out--;
            s->strm->total_out_lo32++;
            if(s->strm->total_out_lo32 == 0) s->strm->total_out_hi32++;
         }
         /*can a new run be started?*/
         if(s->nblock_used == s->save_nblock + 1) 
         return;
         s->state_out_len = 1;
         s->state_out_ch = s->k0;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 2;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 3;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         k1 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
         s->state_out_len = ((Int32) k1) + 4;
         clava_dcg_global[ 63 ]++;
         s->k0 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         if(s->rNToGo == 0) {
            s->rNToGo = BZ2_rNums[s->rTPos];
            s->rTPos++;
            if(s->rTPos == 512) s->rTPos = 0;
         }
         s->rNToGo--;
         ;
         s->k0 ^= ((s->rNToGo == 1) ? 1 : 0);
         s->nblock_used++;
      }
   }
   else {
      while(((Bool) 1)) {
         /*try to finish existing run*/
         while(((Bool) 1)) {
            if(s->strm->avail_out == 0) 
            return;
            if(s->state_out_len == 0) break;
            *((UChar *) (s->strm->next_out)) = s->state_out_ch;
            {
               s->calculatedBlockCRC = (s->calculatedBlockCRC << 8) ^ BZ2_crc32Table[(s->calculatedBlockCRC >> 24) ^ ((UChar) s->state_out_ch)];
            }
            ;
            s->state_out_len--;
            s->strm->next_out++;
            s->strm->avail_out--;
            s->strm->total_out_lo32++;
            if(s->strm->total_out_lo32 == 0) s->strm->total_out_hi32++;
         }
         /*can a new run be started?*/
         if(s->nblock_used == s->save_nblock + 1) 
         return;
         s->state_out_len = 1;
         s->state_out_ch = s->k0;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 2;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         s->state_out_len = 3;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         s->nblock_used++;
         if(s->nblock_used == s->save_nblock + 1) continue;
         if(k1 != s->k0) {
            s->k0 = k1;
            continue;
         }
         ;
         clava_dcg_global[ 63 ]++;
         k1 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         s->nblock_used++;
         s->state_out_len = ((Int32) k1) + 4;
         clava_dcg_global[ 63 ]++;
         s->k0 = BZ2_indexIntoF(s->tPos, s->cftab);
         s->tPos = (((UInt32) s->ll16[s->tPos]) | (((((UInt32) (s->ll4[(s->tPos) >> 1])) >> (((s->tPos) << 2) & 0x4)) & 0xF) << 16));
         ;
         s->nblock_used++;
      }
   }
}

/*---------------------------------------------------*/
int BZ2_bzDecompress(bz_stream *strm) {
   DState *s;
   if(strm == ((void *) 0)) 
   return (-2);
   s = strm->state;
   if(s == ((void *) 0)) 
   return (-2);
   if(s->strm != strm) 
   return (-2);
   while(((Bool) 1)) {
      if(s->state == 1) 
      return (-1);
      if(s->state == 2) {
         if(s->smallDecompress) {
            clava_dcg_global[ 64 ]++;
            unRLE_obuf_to_output_SMALL(s);
         }
         else {
            clava_dcg_global[ 65 ]++;
            unRLE_obuf_to_output_FAST(s);
         }
         if(s->nblock_used == s->save_nblock + 1 && s->state_out_len == 0) {
            {
               s->calculatedBlockCRC = ~(s->calculatedBlockCRC);
            }
            ;
            if(s->verbosity >= 3) {
               clava_dcg_global[ 66 ]++;
               fprintf(stderr, " {0x%x, 0x%x}", s->storedBlockCRC, s->calculatedBlockCRC);
            }
            if(s->verbosity >= 2) {
               clava_dcg_global[ 66 ]++;
               fprintf(stderr, "]");
            }
            if(s->calculatedBlockCRC != s->storedBlockCRC) 
            return (-4);
            s->calculatedCombinedCRC = (s->calculatedCombinedCRC << 1) | (s->calculatedCombinedCRC >> 31);
            s->calculatedCombinedCRC ^= s->calculatedBlockCRC;
            s->state = 14;
         }
         else {
            
            return 0;
         }
      }
      if(s->state >= 10) {
         clava_dcg_global[ 67 ]++;
         Int32 r = BZ2_decompress(s);
         if(r == 4) {
            if(s->verbosity >= 3) {
               clava_dcg_global[ 66 ]++;
               fprintf(stderr, "\n    combined CRCs: stored = 0x%x, computed = 0x%x", s->storedCombinedCRC, s->calculatedCombinedCRC);
            }
            if(s->calculatedCombinedCRC != s->storedCombinedCRC) 
            return (-4);
            
            return r;
         }
         if(s->state != 2) 
         return r;
      }
   }
   {
      if(!(0)) {
         clava_dcg_global[ 68 ]++;
         BZ2_bz__AssertH__fail(6001);
      }
   }
   ;
   
   return 0;
   /*NOTREACHED*/
}

/*---------------------------------------------------*/
int BZ2_bzDecompressEnd(bz_stream *strm) {
   DState *s;
   if(strm == ((void *) 0)) 
   return (-2);
   s = strm->state;
   if(s == ((void *) 0)) 
   return (-2);
   if(s->strm != strm) 
   return (-2);
   if(s->tt != ((void *) 0)) {
      clava_dcg_global[ 69 ]++;
      (strm->bzfree)(strm->opaque, (s->tt));
   }
   if(s->ll16 != ((void *) 0)) {
      clava_dcg_global[ 69 ]++;
      (strm->bzfree)(strm->opaque, (s->ll16));
   }
   if(s->ll4 != ((void *) 0)) {
      clava_dcg_global[ 69 ]++;
      (strm->bzfree)(strm->opaque, (s->ll4));
   }
   clava_dcg_global[ 69 ]++;
   (strm->bzfree)(strm->opaque, (strm->state));
   strm->state = ((void *) 0);
   
   return 0;
}

/*---------------------------------------------------*/
/*--- File I/O stuff                              ---*/
/*---------------------------------------------------*/

struct anon_bzip2_c_4295 {
   FILE *handle;
   Char buf[5000];
   Int32 bufN;
   Bool writing;
   bz_stream strm;
   Int32 lastErr;
   Bool initialisedOk;
};

typedef struct anon_bzip2_c_4295 bzFile;
/*---------------------------------------------*/
static Bool myfeof(FILE *f) {
   clava_dcg_global[ 70 ]++;
   Int32 c = fgetc(f);
   if(c == (-1)) 
   return ((Bool) 1);
   clava_dcg_global[ 71 ]++;
   ungetc(c, f);
   
   return ((Bool) 0);
}

/*---------------------------------------------------*/
BZFILE * BZ2_bzWriteOpen(int *bzerror, FILE *f, int blockSize100k, int verbosity, int workFactor) {
   Int32 ret;
   bzFile *bzf = ((void *) 0);
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   if(f == ((void *) 0) || (blockSize100k < 1 || blockSize100k > 9) || (workFactor < 0 || workFactor > 250) || (verbosity < 0 || verbosity > 4)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   clava_dcg_global[ 72 ]++;
   if(ferror(f)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-6);
         if(bzf != ((void *) 0)) bzf->lastErr = (-6);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   clava_dcg_global[ 73 ]++;
   bzf = malloc(sizeof(bzFile));
   if(bzf == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-3);
         if(bzf != ((void *) 0)) bzf->lastErr = (-3);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   bzf->initialisedOk = ((Bool) 0);
   bzf->bufN = 0;
   bzf->handle = f;
   bzf->writing = ((Bool) 1);
   bzf->strm.bzalloc = ((void *) 0);
   bzf->strm.bzfree = ((void *) 0);
   bzf->strm.opaque = ((void *) 0);
   if(workFactor == 0) workFactor = 30;
   clava_dcg_global[ 74 ]++;
   ret = BZ2_bzCompressInit(&(bzf->strm), blockSize100k, verbosity, workFactor);
   if(ret != 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = ret;
         if(bzf != ((void *) 0)) bzf->lastErr = ret;
      }
      ;
      clava_dcg_global[ 75 ]++;
      free(bzf);
      
      return ((void *) 0);
   }
   ;
   bzf->strm.avail_in = 0;
   bzf->initialisedOk = ((Bool) 1);
   
   return bzf;
}

/*---------------------------------------------------*/
void BZ2_bzWrite(int *bzerror, BZFILE *b, void *buf, int len) {
   Int32 n, n2, ret;
   bzFile *bzf = (bzFile *) b;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   if(bzf == ((void *) 0) || buf == ((void *) 0) || len < 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return;
   }
   ;
   if(!(bzf->writing)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-1);
         if(bzf != ((void *) 0)) bzf->lastErr = (-1);
      }
      ;
      
      return;
   }
   ;
   clava_dcg_global[ 76 ]++;
   if(ferror(bzf->handle)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-6);
         if(bzf != ((void *) 0)) bzf->lastErr = (-6);
      }
      ;
      
      return;
   }
   ;
   if(len == 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = 0;
         if(bzf != ((void *) 0)) bzf->lastErr = 0;
      }
      ;
      
      return;
   }
   ;
   bzf->strm.avail_in = len;
   bzf->strm.next_in = buf;
   while(((Bool) 1)) {
      bzf->strm.avail_out = 5000;
      bzf->strm.next_out = bzf->buf;
      clava_dcg_global[ 77 ]++;
      ret = BZ2_bzCompress(&(bzf->strm), 0);
      if(ret != 1) {
         {
            if(bzerror != ((void *) 0)) *bzerror = ret;
            if(bzf != ((void *) 0)) bzf->lastErr = ret;
         }
         ;
         
         return;
      }
      ;
      if(bzf->strm.avail_out < 5000) {
         n = 5000 - bzf->strm.avail_out;
         clava_dcg_global[ 78 ]++;
         n2 = fwrite((void *) (bzf->buf), sizeof(UChar), n, bzf->handle);
         clava_dcg_global[ 76 ]++;
         if(n != n2 || ferror(bzf->handle)) {
            {
               if(bzerror != ((void *) 0)) *bzerror = (-6);
               if(bzf != ((void *) 0)) bzf->lastErr = (-6);
            }
            ;
            
            return;
         }
         ;
      }
      if(bzf->strm.avail_in == 0) {
         {
            if(bzerror != ((void *) 0)) *bzerror = 0;
            if(bzf != ((void *) 0)) bzf->lastErr = 0;
         }
         ;
         
         return;
      }
      ;
   }
}

/*---------------------------------------------------*/
void BZ2_bzWriteClose(int *bzerror, BZFILE *b, int abandon, unsigned int *nbytes_in, unsigned int *nbytes_out) {
   clava_dcg_global[ 79 ]++;
   BZ2_bzWriteClose64(bzerror, b, abandon, nbytes_in, ((void *) 0), nbytes_out, ((void *) 0));
}

void BZ2_bzWriteClose64(int *bzerror, BZFILE *b, int abandon, unsigned int *nbytes_in_lo32, unsigned int *nbytes_in_hi32, unsigned int *nbytes_out_lo32, unsigned int *nbytes_out_hi32) {
   Int32 n, n2, ret;
   bzFile *bzf = (bzFile *) b;
   if(bzf == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = 0;
         if(bzf != ((void *) 0)) bzf->lastErr = 0;
      }
      ;
      
      return;
   }
   ;
   if(!(bzf->writing)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-1);
         if(bzf != ((void *) 0)) bzf->lastErr = (-1);
      }
      ;
      
      return;
   }
   ;
   clava_dcg_global[ 80 ]++;
   if(ferror(bzf->handle)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-6);
         if(bzf != ((void *) 0)) bzf->lastErr = (-6);
      }
      ;
      
      return;
   }
   ;
   if(nbytes_in_lo32 != ((void *) 0)) *nbytes_in_lo32 = 0;
   if(nbytes_in_hi32 != ((void *) 0)) *nbytes_in_hi32 = 0;
   if(nbytes_out_lo32 != ((void *) 0)) *nbytes_out_lo32 = 0;
   if(nbytes_out_hi32 != ((void *) 0)) *nbytes_out_hi32 = 0;
   if((!abandon) && bzf->lastErr == 0) {
      while(((Bool) 1)) {
         bzf->strm.avail_out = 5000;
         bzf->strm.next_out = bzf->buf;
         clava_dcg_global[ 81 ]++;
         ret = BZ2_bzCompress(&(bzf->strm), 2);
         if(ret != 3 && ret != 4) {
            {
               if(bzerror != ((void *) 0)) *bzerror = ret;
               if(bzf != ((void *) 0)) bzf->lastErr = ret;
            }
            ;
            
            return;
         }
         ;
         if(bzf->strm.avail_out < 5000) {
            n = 5000 - bzf->strm.avail_out;
            clava_dcg_global[ 82 ]++;
            n2 = fwrite((void *) (bzf->buf), sizeof(UChar), n, bzf->handle);
            clava_dcg_global[ 80 ]++;
            if(n != n2 || ferror(bzf->handle)) {
               {
                  if(bzerror != ((void *) 0)) *bzerror = (-6);
                  if(bzf != ((void *) 0)) bzf->lastErr = (-6);
               }
               ;
               
               return;
            }
            ;
         }
         if(ret == 4) break;
      }
   }
   clava_dcg_global[ 80 ]++;
   if(!abandon && !ferror(bzf->handle)) {
      clava_dcg_global[ 83 ]++;
      fflush(bzf->handle);
      clava_dcg_global[ 80 ]++;
      if(ferror(bzf->handle)) {
         {
            if(bzerror != ((void *) 0)) *bzerror = (-6);
            if(bzf != ((void *) 0)) bzf->lastErr = (-6);
         }
         ;
         
         return;
      }
      ;
   }
   if(nbytes_in_lo32 != ((void *) 0)) *nbytes_in_lo32 = bzf->strm.total_in_lo32;
   if(nbytes_in_hi32 != ((void *) 0)) *nbytes_in_hi32 = bzf->strm.total_in_hi32;
   if(nbytes_out_lo32 != ((void *) 0)) *nbytes_out_lo32 = bzf->strm.total_out_lo32;
   if(nbytes_out_hi32 != ((void *) 0)) *nbytes_out_hi32 = bzf->strm.total_out_hi32;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   clava_dcg_global[ 84 ]++;
   BZ2_bzCompressEnd(&(bzf->strm));
   clava_dcg_global[ 85 ]++;
   free(bzf);
}

/*---------------------------------------------------*/
BZFILE * BZ2_bzReadOpen(int *bzerror, FILE *f, int verbosity, int small, void *unused, int nUnused) {
   bzFile *bzf = ((void *) 0);
   int ret;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   if(f == ((void *) 0) || (small != 0 && small != 1) || (verbosity < 0 || verbosity > 4) || (unused == ((void *) 0) && nUnused != 0) || (unused != ((void *) 0) && (nUnused < 0 || nUnused > 5000))) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   clava_dcg_global[ 86 ]++;
   if(ferror(f)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-6);
         if(bzf != ((void *) 0)) bzf->lastErr = (-6);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   clava_dcg_global[ 87 ]++;
   bzf = malloc(sizeof(bzFile));
   if(bzf == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-3);
         if(bzf != ((void *) 0)) bzf->lastErr = (-3);
      }
      ;
      
      return ((void *) 0);
   }
   ;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   bzf->initialisedOk = ((Bool) 0);
   bzf->handle = f;
   bzf->bufN = 0;
   bzf->writing = ((Bool) 0);
   bzf->strm.bzalloc = ((void *) 0);
   bzf->strm.bzfree = ((void *) 0);
   bzf->strm.opaque = ((void *) 0);
   while(nUnused > 0) {
      bzf->buf[bzf->bufN] = *((UChar *) (unused));
      bzf->bufN++;
      unused = ((void *) (1 + ((UChar *) (unused))));
      nUnused--;
   }
   clava_dcg_global[ 88 ]++;
   ret = BZ2_bzDecompressInit(&(bzf->strm), verbosity, small);
   if(ret != 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = ret;
         if(bzf != ((void *) 0)) bzf->lastErr = ret;
      }
      ;
      clava_dcg_global[ 89 ]++;
      free(bzf);
      
      return ((void *) 0);
   }
   ;
   bzf->strm.avail_in = bzf->bufN;
   bzf->strm.next_in = bzf->buf;
   bzf->initialisedOk = ((Bool) 1);
   
   return bzf;
}

/*---------------------------------------------------*/
void BZ2_bzReadClose(int *bzerror, BZFILE *b) {
   bzFile *bzf = (bzFile *) b;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   if(bzf == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = 0;
         if(bzf != ((void *) 0)) bzf->lastErr = 0;
      }
      ;
      
      return;
   }
   ;
   if(bzf->writing) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-1);
         if(bzf != ((void *) 0)) bzf->lastErr = (-1);
      }
      ;
      
      return;
   }
   ;
   if(bzf->initialisedOk) {
      clava_dcg_global[ 90 ]++;
      (void) BZ2_bzDecompressEnd(&(bzf->strm));
   }
   clava_dcg_global[ 91 ]++;
   free(bzf);
}

/*---------------------------------------------------*/
int BZ2_bzRead(int *bzerror, BZFILE *b, void *buf, int len) {
   Int32 n, ret;
   bzFile *bzf = (bzFile *) b;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   if(bzf == ((void *) 0) || buf == ((void *) 0) || len < 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return 0;
   }
   ;
   if(bzf->writing) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-1);
         if(bzf != ((void *) 0)) bzf->lastErr = (-1);
      }
      ;
      
      return 0;
   }
   ;
   if(len == 0) {
      {
         if(bzerror != ((void *) 0)) *bzerror = 0;
         if(bzf != ((void *) 0)) bzf->lastErr = 0;
      }
      ;
      
      return 0;
   }
   ;
   bzf->strm.avail_out = len;
   bzf->strm.next_out = buf;
   while(((Bool) 1)) {
      clava_dcg_global[ 92 ]++;
      if(ferror(bzf->handle)) {
         {
            if(bzerror != ((void *) 0)) *bzerror = (-6);
            if(bzf != ((void *) 0)) bzf->lastErr = (-6);
         }
         ;
         
         return 0;
      }
      ;
      clava_dcg_global[ 93 ]++;
      if(bzf->strm.avail_in == 0 && !myfeof(bzf->handle)) {
         clava_dcg_global[ 94 ]++;
         n = fread(bzf->buf, sizeof(UChar), 5000, bzf->handle);
         clava_dcg_global[ 92 ]++;
         if(ferror(bzf->handle)) {
            {
               if(bzerror != ((void *) 0)) *bzerror = (-6);
               if(bzf != ((void *) 0)) bzf->lastErr = (-6);
            }
            ;
            
            return 0;
         }
         ;
         bzf->bufN = n;
         bzf->strm.avail_in = bzf->bufN;
         bzf->strm.next_in = bzf->buf;
      }
      clava_dcg_global[ 95 ]++;
      ret = BZ2_bzDecompress(&(bzf->strm));
      if(ret != 0 && ret != 4) {
         {
            if(bzerror != ((void *) 0)) *bzerror = ret;
            if(bzf != ((void *) 0)) bzf->lastErr = ret;
         }
         ;
         
         return 0;
      }
      ;
      clava_dcg_global[ 93 ]++;
      if(ret == 0 && myfeof(bzf->handle) && bzf->strm.avail_in == 0 && bzf->strm.avail_out > 0) {
         {
            if(bzerror != ((void *) 0)) *bzerror = (-7);
            if(bzf != ((void *) 0)) bzf->lastErr = (-7);
         }
         ;
         
         return 0;
      }
      ;
      if(ret == 4) {
         {
            if(bzerror != ((void *) 0)) *bzerror = 4;
            if(bzf != ((void *) 0)) bzf->lastErr = 4;
         }
         ;
         
         return len - bzf->strm.avail_out;
      }
      ;
      if(bzf->strm.avail_out == 0) {
         {
            if(bzerror != ((void *) 0)) *bzerror = 0;
            if(bzf != ((void *) 0)) bzf->lastErr = 0;
         }
         ;
         
         return len;
      }
      ;
   }
   
   return 0;
   /*not reached*/
}

/*---------------------------------------------------*/
void BZ2_bzReadGetUnused(int *bzerror, BZFILE *b, void **unused, int *nUnused) {
   bzFile *bzf = (bzFile *) b;
   if(bzf == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return;
   }
   ;
   if(bzf->lastErr != 4) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-1);
         if(bzf != ((void *) 0)) bzf->lastErr = (-1);
      }
      ;
      
      return;
   }
   ;
   if(unused == ((void *) 0) || nUnused == ((void *) 0)) {
      {
         if(bzerror != ((void *) 0)) *bzerror = (-2);
         if(bzf != ((void *) 0)) bzf->lastErr = (-2);
      }
      ;
      
      return;
   }
   ;
   {
      if(bzerror != ((void *) 0)) *bzerror = 0;
      if(bzf != ((void *) 0)) bzf->lastErr = 0;
   }
   ;
   *nUnused = bzf->strm.avail_in;
   *unused = bzf->strm.next_in;
}

/*---------------------------------------------------*/
/*--- Misc convenience stuff                      ---*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
int BZ2_bzBuffToBuffCompress(char *dest, unsigned int *destLen, char *source, unsigned int sourceLen, int blockSize100k, int verbosity, int workFactor) {
   bz_stream strm;
   int ret;
   if(dest == ((void *) 0) || destLen == ((void *) 0) || source == ((void *) 0) || blockSize100k < 1 || blockSize100k > 9 || verbosity < 0 || verbosity > 4 || workFactor < 0 || workFactor > 250) 
   return (-2);
   if(workFactor == 0) workFactor = 30;
   strm.bzalloc = ((void *) 0);
   strm.bzfree = ((void *) 0);
   strm.opaque = ((void *) 0);
   clava_dcg_global[ 96 ]++;
   ret = BZ2_bzCompressInit(&strm, blockSize100k, verbosity, workFactor);
   if(ret != 0) 
   return ret;
   strm.next_in = source;
   strm.next_out = dest;
   strm.avail_in = sourceLen;
   strm.avail_out = *destLen;
   clava_dcg_global[ 97 ]++;
   ret = BZ2_bzCompress(&strm, 2);
   if(ret == 3) goto output_overflow;
   if(ret != 4) goto errhandler;
   /*normal termination*/
   *destLen -= strm.avail_out;
   clava_dcg_global[ 98 ]++;
   BZ2_bzCompressEnd(&strm);
   
   return 0;
   output_overflow:
   clava_dcg_global[ 98 ]++;
   BZ2_bzCompressEnd(&strm);
   
   return (-8);
   errhandler:
   clava_dcg_global[ 98 ]++;
   BZ2_bzCompressEnd(&strm);
   
   return ret;
}

/*---------------------------------------------------*/
int BZ2_bzBuffToBuffDecompress(char *dest, unsigned int *destLen, char *source, unsigned int sourceLen, int small, int verbosity) {
   bz_stream strm;
   int ret;
   if(dest == ((void *) 0) || destLen == ((void *) 0) || source == ((void *) 0) || (small != 0 && small != 1) || verbosity < 0 || verbosity > 4) 
   return (-2);
   strm.bzalloc = ((void *) 0);
   strm.bzfree = ((void *) 0);
   strm.opaque = ((void *) 0);
   clava_dcg_global[ 99 ]++;
   ret = BZ2_bzDecompressInit(&strm, verbosity, small);
   if(ret != 0) 
   return ret;
   strm.next_in = source;
   strm.next_out = dest;
   strm.avail_in = sourceLen;
   strm.avail_out = *destLen;
   clava_dcg_global[ 100 ]++;
   ret = BZ2_bzDecompress(&strm);
   if(ret == 0) goto output_overflow_or_eof;
   if(ret != 4) goto errhandler;
   /*normal termination*/
   *destLen -= strm.avail_out;
   clava_dcg_global[ 101 ]++;
   BZ2_bzDecompressEnd(&strm);
   
   return 0;
   output_overflow_or_eof:
   if(strm.avail_out > 0) {
      clava_dcg_global[ 101 ]++;
      BZ2_bzDecompressEnd(&strm);
      
      return (-7);
   }
   else {
      clava_dcg_global[ 101 ]++;
      BZ2_bzDecompressEnd(&strm);
      
      return (-8);
   }
   ;
   errhandler:
   clava_dcg_global[ 101 ]++;
   BZ2_bzDecompressEnd(&strm);
   
   return ret;
}

/*---------------------------------------------------*/
/*--
Code contributed by Yoshioka Tsuneo
(QWF00133@niftyserve.or.jp/tsuneo-y@is.aist-nara.ac.jp),
to support better zlib compatibility.
This code is not _officially_ part of libbzip2 (yet);
I haven't tested it, documented it, or considered the
threading-safeness of it.
If this code breaks, please contact both Yoshioka and me.
--*/
/*---------------------------------------------------*/
/*---------------------------------------------------*/
/*--
return version like "0.9.0c".
--*/
char const * BZ2_bzlibVersion() {
   
   return "1.0.2, 30-Dec-2001";
}

/*---------------------------------------------------*/
static BZFILE * bzopen_or_bzdopen(char const *path, int fd, char const *mode, int open_mode) {
   /*no use when bzdopen*/
   /*no use when bzdopen*/
   /*bzopen: 0, bzdopen:1*/
   int bzerr;
   char unused[5000];
   int blockSize100k = 9;
   int writing = 0;
   char mode2[10] = "";
   FILE *fp = ((void *) 0);
   BZFILE *bzfp = ((void *) 0);
   int verbosity = 0;
   int workFactor = 30;
   int smallMode = 0;
   int nUnused = 0;
   if(mode == ((void *) 0)) 
   return ((void *) 0);
   while(*mode) {
      switch (*mode) {
         case 'r':
         writing = 0;
         break;
         case 'w':
         writing = 1;
         break;
         case 's':
         smallMode = 1;
         break;
         default:
         clava_dcg_global[ 102 ]++;
         if(((*__ctype_b_loc())[(int) (((int) (*mode)))] & (unsigned short) _ISdigit)) {
            blockSize100k = *mode - 0x30;
         }
      }
      mode++;
   }
   clava_dcg_global[ 103 ]++;
   strcat(mode2, writing ? "w" : "r");
   clava_dcg_global[ 103 ]++;
   strcat(mode2, "b");
   if(open_mode == 0) {
      clava_dcg_global[ 104 ]++;
      if(path == ((void *) 0) || strcmp(path, "") == 0) {
         fp = (writing ? stdout : stdin);
      }
      else {
         clava_dcg_global[ 105 ]++;
         fp = fopen(path, mode2);
      }
   }
   else {
      clava_dcg_global[ 106 ]++;
      fp = fdopen(fd, mode2);
   }
   if(fp == ((void *) 0)) 
   return ((void *) 0);
   if(writing) {
      /*Guard against total chaos and anarchy -- JRS*/
      if(blockSize100k < 1) blockSize100k = 1;
      if(blockSize100k > 9) blockSize100k = 9;
      clava_dcg_global[ 107 ]++;
      bzfp = BZ2_bzWriteOpen(&bzerr, fp, blockSize100k, verbosity, workFactor);
   }
   else {
      clava_dcg_global[ 108 ]++;
      bzfp = BZ2_bzReadOpen(&bzerr, fp, verbosity, smallMode, unused, nUnused);
   }
   if(bzfp == ((void *) 0)) {
      if(fp != stdin && fp != stdout) {
         clava_dcg_global[ 109 ]++;
         fclose(fp);
      }
      
      return ((void *) 0);
   }
   
   return bzfp;
}

/*---------------------------------------------------*/
/*--
open file for read or write.
ex) bzopen("file","w9")
case path="" or NULL => use stdin or stdout.
--*/
BZFILE * BZ2_bzopen(char const *path, char const *mode) {
   clava_dcg_global[ 110 ]++;
   
   return bzopen_or_bzdopen(path, -1, mode, 0);
   /*bzopen*/
}

/*---------------------------------------------------*/
BZFILE * BZ2_bzdopen(int fd, char const *mode) {
   clava_dcg_global[ 111 ]++;
   
   return bzopen_or_bzdopen(((void *) 0), fd, mode, 1);
   /*bzdopen*/
}

/*---------------------------------------------------*/
int BZ2_bzread(BZFILE *b, void *buf, int len) {
   int bzerr, nread;
   if(((bzFile *) b)->lastErr == 4) 
   return 0;
   clava_dcg_global[ 112 ]++;
   nread = BZ2_bzRead(&bzerr, b, buf, len);
   if(bzerr == 0 || bzerr == 4) {
      
      return nread;
   }
   else {
      
      return -1;
   }
}

/*---------------------------------------------------*/
int BZ2_bzwrite(BZFILE *b, void *buf, int len) {
   int bzerr;
   clava_dcg_global[ 113 ]++;
   BZ2_bzWrite(&bzerr, b, buf, len);
   if(bzerr == 0) {
      
      return len;
   }
   else {
      
      return -1;
   }
}

/*---------------------------------------------------*/
int BZ2_bzflush(BZFILE *b) {
   /*do nothing now...*/
   
   return 0;
}

/*---------------------------------------------------*/
void BZ2_bzclose(BZFILE *b) {
   int bzerr;
   FILE *fp = ((bzFile *) b)->handle;
   if(b == ((void *) 0)) {
      
      return;
   }
   if(((bzFile *) b)->writing) {
      clava_dcg_global[ 114 ]++;
      BZ2_bzWriteClose(&bzerr, b, 0, ((void *) 0), ((void *) 0));
      if(bzerr != 0) {
         clava_dcg_global[ 114 ]++;
         BZ2_bzWriteClose(((void *) 0), b, 1, ((void *) 0), ((void *) 0));
      }
   }
   else {
      clava_dcg_global[ 115 ]++;
      BZ2_bzReadClose(&bzerr, b);
   }
   if(fp != stdin && fp != stdout) {
      clava_dcg_global[ 116 ]++;
      fclose(fp);
   }
}

/*---------------------------------------------------*/

/*--
return last error code
--*/

static char *bzerrorstrings[16] = {"OK", "SEQUENCE_ERROR", "PARAM_ERROR", "MEM_ERROR", "DATA_ERROR", "DATA_ERROR_MAGIC", "IO_ERROR", "UNEXPECTED_EOF", "OUTBUFF_FULL", "CONFIG_ERROR", "???", "???", "???", "???", "???", "???"};
char const * BZ2_bzerror(BZFILE *b, int *errnum) {
   int err = ((bzFile *) b)->lastErr;
   if(err > 0) err = 0;
   *errnum = err;
   
   return bzerrorstrings[err * -1];
}

/*-------------------------------------------------------------*/

/*--- end                                           bzlib.c ---*/

/*-------------------------------------------------------------*/

/*-----------------------------------------------------------*/

/*--- A block-sorting, lossless compressor        bzip2.c ---*/

/*-----------------------------------------------------------*/

/*----------------------------------------------------*/

/*--- IMPORTANT                                    ---*/

/*----------------------------------------------------*/

/*--
WARNING:
This program and library (attempts to) compress data by
performing several non-trivial transformations on it.
Unless you are 100% familiar with *all* the algorithms
contained herein, and with the consequences of modifying them,
you should NOT meddle with the compression or decompression
machinery.  Incorrect changes can and very likely *will*
lead to disasterous loss of data.

DISCLAIMER:
I TAKE NO RESPONSIBILITY FOR ANY LOSS OF DATA ARISING FROM THE
USE OF THIS PROGRAM, HOWSOEVER CAUSED.

Every compression of a file implies an assumption that the
compressed file can be decompressed to reproduce the original.
Great efforts in design, coding and testing have been made to
ensure that this program works correctly.  However, the
complexity of the algorithms, and, in particular, the presence
of various special cases in the code which occur with very low
but non-zero probability make it impossible to rule out the
possibility of bugs remaining in the program.  DO NOT COMPRESS
ANY DATA WITH THIS PROGRAM AND/OR LIBRARY UNLESS YOU ARE PREPARED
TO ACCEPT THE POSSIBILITY, HOWEVER SMALL, THAT THE DATA WILL
NOT BE RECOVERABLE.

That is not to say this program is inherently unreliable.
Indeed, I very much hope the opposite is true.  bzip2/libbzip2
has been carefully constructed and extensively tested.

PATENTS:
To the best of my knowledge, bzip2/libbzip2 does not use any
patented algorithms.  However, I do not have the resources
available to carry out a full patent search.  Therefore I cannot
give any guarantee of the above statement.
--*/

/*----------------------------------------------------*/

/*--- and now for something much more pleasant :-) ---*/

/*----------------------------------------------------*/

/*---------------------------------------------*/

/*--
Place a 1 beside your platform, and 0 elsewhere.
--*/

/*--
Generic 32-bit Unix.
Also works on 64-bit Unix boxes.
This is the default.
--*/

/*--
Win32, as seen by Jacob Navia's excellent
port of (Chris Fraser & David Hanson)'s excellent
lcc compiler.  Or with MS Visual C.
This is selected automatically if compiled by a compiler which
defines _WIN32, not including the Cygwin GCC.
--*/

/*---------------------------------------------*/

/*--
Some stuff for all platforms.
--*/

/*---------------------------------------------*/

/*--
Platform-specific stuff.
--*/

/**/

/**/

/*BZ_UNIX*/

/**/

/*BZ_LCCWIN32*/

/*---------------------------------------------*/

/*--
Some more stuff for all platforms :-)
--*/

/*--
IntNative is your platform's `native' int size.
Only here to avoid probs with 64-bit platforms.
--*/

typedef int IntNative;
/*---------------------------------------------------*/

/*--- Misc (file handling) data decls             ---*/

/*---------------------------------------------------*/

Int32 verbosity;
Bool keepInputFiles;
Bool smallMode;
Bool deleteOutputOnInterrupt;
Bool forceOverwrite;
Bool testFailsExist;
Bool unzFailsExist;
Bool noisy;
Int32 numFileNames;
Int32 numFilesProcessed;
Int32 blockSize100k;
Int32 exitValue;
/*-- source modes; F==file, I==stdin, O==stdout --*/

/*-- operation modes --*/

Int32 opMode;
Int32 srcMode;
Int32 longestFileName;
Char inName[1034];
Char outName[1034];
Char tmpName[1034];
Char *progName;
Char progNameReally[1034];
FILE *outputHandleJustInCase;
Int32 workFactor;
static void panic(Char *);
static void ioError();
static void outOfMemory();
static void configError();
static void crcError();
static void cleanUpAndFail(Int32);
static void compressedStreamEOF();
static void copyFileName(Char *, Char *);
static void * myMalloc(Int32);
/*---------------------------------------------------*/
/*--- An implementation of 64-bit ints.  Sigh.    ---*/
/*--- Roll on widespread deployment of ANSI C9X ! ---*/
/*---------------------------------------------------*/

struct anon_bzip2_c_5228 {
   UChar b[8];
};

typedef struct anon_bzip2_c_5228 UInt64;
static void uInt64_from_UInt32s(UInt64 *n, UInt32 lo32, UInt32 hi32) {
   n->b[7] = (UChar) ((hi32 >> 24) & 0xFF);
   n->b[6] = (UChar) ((hi32 >> 16) & 0xFF);
   n->b[5] = (UChar) ((hi32 >> 8) & 0xFF);
   n->b[4] = (UChar) (hi32 & 0xFF);
   n->b[3] = (UChar) ((lo32 >> 24) & 0xFF);
   n->b[2] = (UChar) ((lo32 >> 16) & 0xFF);
   n->b[1] = (UChar) ((lo32 >> 8) & 0xFF);
   n->b[0] = (UChar) (lo32 & 0xFF);
}

static double uInt64_to_double(UInt64 *n) {
   Int32 i;
   double base = 1.0;
   double sum = 0.0;
   for(i = 0; i < 8; i++) {
      sum += base * (double) (n->b[i]);
      base *= 256.0;
   }
   
   return sum;
}

static Bool uInt64_isZero(UInt64 *n) {
   Int32 i;
   for(i = 0; i < 8; i++) if(n->b[i] != 0) 
   return 0;
   
   return 1;
}

/*Divide *n by 10, and return the remainder.*/
static Int32 uInt64_qrm10(UInt64 *n) {
   UInt32 rem, tmp;
   Int32 i;
   rem = 0;
   for(i = 7; i >= 0; i--) {
      tmp = rem * 256 + n->b[i];
      n->b[i] = tmp / 10;
      rem = tmp % 10;
   }
   
   return rem;
}

/*... and the Whole Entire Point of all this UInt64 stuff is
so that we can supply the following function.
*/
static void uInt64_toAscii(char *outbuf, UInt64 *n) {
   Int32 i, q;
   UChar buf[32];
   Int32 nBuf = 0;
   UInt64 n_copy = *n;
   do  {
      clava_dcg_global[ 118 ]++;
      clava_dcg_global[ 117 ]++;
      q = uInt64_qrm10(&n_copy);
      buf[nBuf] = q + '0';
      nBuf++;
   }
   while (!uInt64_isZero(&n_copy));
   outbuf[nBuf] = 0;
   for(i = 0; i < nBuf; i++) outbuf[i] = buf[nBuf - i - 1];
}

/*---------------------------------------------------*/
/*--- Processing of complete files and streams    ---*/
/*---------------------------------------------------*/
/*---------------------------------------------*/
/*---------------------------------------------*/
static void compressStream(FILE *stream, FILE *zStream) {
   BZFILE *bzf = ((void *) 0);
   UChar ibuf[5000];
   Int32 nIbuf;
   UInt32 nbytes_in_lo32, nbytes_in_hi32;
   UInt32 nbytes_out_lo32, nbytes_out_hi32;
   Int32 bzerr, bzerr_dummy, ret;
   ;
   ;
   clava_dcg_global[ 119 ]++;
   if(ferror(stream)) goto errhandler_io;
   clava_dcg_global[ 119 ]++;
   if(ferror(zStream)) goto errhandler_io;
   clava_dcg_global[ 120 ]++;
   bzf = BZ2_bzWriteOpen(&bzerr, zStream, blockSize100k, verbosity, workFactor);
   if(bzerr != 0) goto errhandler;
   if(verbosity >= 2) {
      clava_dcg_global[ 121 ]++;
      fprintf(stderr, "\n");
   }
   while(((Bool) 1)) {
      clava_dcg_global[ 122 ]++;
      if(myfeof(stream)) break;
      clava_dcg_global[ 123 ]++;
      nIbuf = fread(ibuf, sizeof(UChar), 5000, stream);
      clava_dcg_global[ 119 ]++;
      if(ferror(stream)) goto errhandler_io;
      if(nIbuf > 0) {
         clava_dcg_global[ 124 ]++;
         BZ2_bzWrite(&bzerr, bzf, (void *) ibuf, nIbuf);
      }
      if(bzerr != 0) goto errhandler;
   }
   clava_dcg_global[ 125 ]++;
   BZ2_bzWriteClose64(&bzerr, bzf, 0, &nbytes_in_lo32, &nbytes_in_hi32, &nbytes_out_lo32, &nbytes_out_hi32);
   if(bzerr != 0) goto errhandler;
   clava_dcg_global[ 119 ]++;
   if(ferror(zStream)) goto errhandler_io;
   clava_dcg_global[ 126 ]++;
   ret = fflush(zStream);
   if(ret == (-1)) goto errhandler_io;
   if(zStream != stdout) {
      clava_dcg_global[ 127 ]++;
      ret = fclose(zStream);
      outputHandleJustInCase = ((void *) 0);
      if(ret == (-1)) goto errhandler_io;
   }
   outputHandleJustInCase = ((void *) 0);
   clava_dcg_global[ 119 ]++;
   if(ferror(stream)) goto errhandler_io;
   clava_dcg_global[ 127 ]++;
   ret = fclose(stream);
   if(ret == (-1)) goto errhandler_io;
   if(verbosity >= 1) {
      if(nbytes_in_lo32 == 0 && nbytes_in_hi32 == 0) {
         clava_dcg_global[ 121 ]++;
         fprintf(stderr, " no data compressed.\n");
      }
      else {
         Char buf_nin[32];
         Char buf_nout[32];
         UInt64 nbytes_in, nbytes_out;
         double nbytes_in_d, nbytes_out_d;
         clava_dcg_global[ 128 ]++;
         uInt64_from_UInt32s(&nbytes_in, nbytes_in_lo32, nbytes_in_hi32);
         clava_dcg_global[ 128 ]++;
         uInt64_from_UInt32s(&nbytes_out, nbytes_out_lo32, nbytes_out_hi32);
         clava_dcg_global[ 129 ]++;
         nbytes_in_d = uInt64_to_double(&nbytes_in);
         clava_dcg_global[ 129 ]++;
         nbytes_out_d = uInt64_to_double(&nbytes_out);
         clava_dcg_global[ 130 ]++;
         uInt64_toAscii(buf_nin, &nbytes_in);
         clava_dcg_global[ 130 ]++;
         uInt64_toAscii(buf_nout, &nbytes_out);
         clava_dcg_global[ 121 ]++;
         fprintf(stderr, "%6.3f:1, %6.3f bits/byte, %5.2f%% saved, %s in, %s out.\n", nbytes_in_d / nbytes_out_d, (8.0 * nbytes_out_d) / nbytes_in_d, 100.0 * (1.0 - nbytes_out_d / nbytes_in_d), buf_nin, buf_nout);
      }
   }
   
   return;
   errhandler:
   clava_dcg_global[ 125 ]++;
   BZ2_bzWriteClose64(&bzerr_dummy, bzf, 1, &nbytes_in_lo32, &nbytes_in_hi32, &nbytes_out_lo32, &nbytes_out_hi32);
   switch (bzerr) {
      case (-9):
      clava_dcg_global[ 131 ]++;
      configError();
      break;
      case (-3):
      clava_dcg_global[ 132 ]++;
      outOfMemory();
      break;
      case (-6):
      errhandler_io:
      clava_dcg_global[ 133 ]++;
      ioError();
      break;
      default:
      clava_dcg_global[ 134 ]++;
      panic("compress:unexpected error");
   }
   clava_dcg_global[ 134 ]++;
   panic("compress:end");
   /*notreached*/
}

/*---------------------------------------------*/
static Bool uncompressStream(FILE *zStream, FILE *stream) {
   BZFILE *bzf = ((void *) 0);
   Int32 bzerr, bzerr_dummy, ret, nread, streamNo, i;
   UChar obuf[5000];
   UChar unused[5000];
   Int32 nUnused;
   UChar *unusedTmp;
   nUnused = 0;
   streamNo = 0;
   ;
   ;
   clava_dcg_global[ 135 ]++;
   if(ferror(stream)) goto errhandler_io;
   clava_dcg_global[ 135 ]++;
   if(ferror(zStream)) goto errhandler_io;
   while(((Bool) 1)) {
      clava_dcg_global[ 136 ]++;
      bzf = BZ2_bzReadOpen(&bzerr, zStream, verbosity, (int) smallMode, unused, nUnused);
      if(bzf == ((void *) 0) || bzerr != 0) goto errhandler;
      streamNo++;
      while(bzerr == 0) {
         clava_dcg_global[ 137 ]++;
         nread = BZ2_bzRead(&bzerr, bzf, obuf, 5000);
         if(bzerr == (-5)) goto trycat;
         if((bzerr == 0 || bzerr == 4) && nread > 0) {
            clava_dcg_global[ 138 ]++;
            fwrite(obuf, sizeof(UChar), nread, stream);
         }
         clava_dcg_global[ 135 ]++;
         if(ferror(stream)) goto errhandler_io;
      }
      if(bzerr != 4) goto errhandler;
      clava_dcg_global[ 139 ]++;
      BZ2_bzReadGetUnused(&bzerr, bzf, (void **) (&unusedTmp), &nUnused);
      if(bzerr != 0) {
         clava_dcg_global[ 140 ]++;
         panic("decompress:bzReadGetUnused");
      }
      for(i = 0; i < nUnused; i++) unused[i] = unusedTmp[i];
      clava_dcg_global[ 141 ]++;
      BZ2_bzReadClose(&bzerr, bzf);
      if(bzerr != 0) {
         clava_dcg_global[ 140 ]++;
         panic("decompress:bzReadGetUnused");
      }
      clava_dcg_global[ 142 ]++;
      if(nUnused == 0 && myfeof(zStream)) break;
   }
   closeok:
   clava_dcg_global[ 135 ]++;
   if(ferror(zStream)) goto errhandler_io;
   clava_dcg_global[ 143 ]++;
   ret = fclose(zStream);
   if(ret == (-1)) goto errhandler_io;
   clava_dcg_global[ 135 ]++;
   if(ferror(stream)) goto errhandler_io;
   clava_dcg_global[ 144 ]++;
   ret = fflush(stream);
   if(ret != 0) goto errhandler_io;
   if(stream != stdout) {
      clava_dcg_global[ 143 ]++;
      ret = fclose(stream);
      outputHandleJustInCase = ((void *) 0);
      if(ret == (-1)) goto errhandler_io;
   }
   outputHandleJustInCase = ((void *) 0);
   if(verbosity >= 2) {
      clava_dcg_global[ 145 ]++;
      fprintf(stderr, "\n    ");
   }
   
   return ((Bool) 1);
   trycat:
   if(forceOverwrite) {
      clava_dcg_global[ 146 ]++;
      rewind(zStream);
      while(((Bool) 1)) {
         clava_dcg_global[ 142 ]++;
         if(myfeof(zStream)) break;
         clava_dcg_global[ 147 ]++;
         nread = fread(obuf, sizeof(UChar), 5000, zStream);
         clava_dcg_global[ 135 ]++;
         if(ferror(zStream)) goto errhandler_io;
         if(nread > 0) {
            clava_dcg_global[ 138 ]++;
            fwrite(obuf, sizeof(UChar), nread, stream);
         }
         clava_dcg_global[ 135 ]++;
         if(ferror(stream)) goto errhandler_io;
      }
      goto closeok;
   }
   errhandler:
   clava_dcg_global[ 141 ]++;
   BZ2_bzReadClose(&bzerr_dummy, bzf);
   switch (bzerr) {
      case (-9):
      clava_dcg_global[ 148 ]++;
      configError();
      break;
      case (-6):
      errhandler_io:
      clava_dcg_global[ 149 ]++;
      ioError();
      break;
      case (-4):
      clava_dcg_global[ 150 ]++;
      crcError();
      case (-3):
      clava_dcg_global[ 151 ]++;
      outOfMemory();
      case (-7):
      clava_dcg_global[ 152 ]++;
      compressedStreamEOF();
      case (-5):
      if(zStream != stdin) {
         clava_dcg_global[ 143 ]++;
         fclose(zStream);
      }
      if(stream != stdout) {
         clava_dcg_global[ 143 ]++;
         fclose(stream);
      }
      if(streamNo == 1) {
         
         return ((Bool) 0);
      }
      else {
         if(noisy) {
            clava_dcg_global[ 145 ]++;
            fprintf(stderr, "\n%s: %s: trailing garbage after EOF ignored\n", progName, inName);
         }
         
         return ((Bool) 1);
      }
      default:
      clava_dcg_global[ 140 ]++;
      panic("decompress:unexpected error");
   }
   clava_dcg_global[ 140 ]++;
   panic("decompress:end");
   
   return ((Bool) 1);
   /*notreached*/
}

/*---------------------------------------------*/
static Bool testStream(FILE *zStream) {
   BZFILE *bzf = ((void *) 0);
   Int32 bzerr, bzerr_dummy, ret, nread, streamNo, i;
   UChar obuf[5000];
   UChar unused[5000];
   Int32 nUnused;
   UChar *unusedTmp;
   nUnused = 0;
   streamNo = 0;
   ;
   clava_dcg_global[ 153 ]++;
   if(ferror(zStream)) goto errhandler_io;
   while(((Bool) 1)) {
      clava_dcg_global[ 154 ]++;
      bzf = BZ2_bzReadOpen(&bzerr, zStream, verbosity, (int) smallMode, unused, nUnused);
      if(bzf == ((void *) 0) || bzerr != 0) goto errhandler;
      streamNo++;
      while(bzerr == 0) {
         clava_dcg_global[ 155 ]++;
         nread = BZ2_bzRead(&bzerr, bzf, obuf, 5000);
         if(bzerr == (-5)) goto errhandler;
      }
      if(bzerr != 4) goto errhandler;
      clava_dcg_global[ 156 ]++;
      BZ2_bzReadGetUnused(&bzerr, bzf, (void **) (&unusedTmp), &nUnused);
      if(bzerr != 0) {
         clava_dcg_global[ 157 ]++;
         panic("test:bzReadGetUnused");
      }
      for(i = 0; i < nUnused; i++) unused[i] = unusedTmp[i];
      clava_dcg_global[ 158 ]++;
      BZ2_bzReadClose(&bzerr, bzf);
      if(bzerr != 0) {
         clava_dcg_global[ 157 ]++;
         panic("test:bzReadGetUnused");
      }
      clava_dcg_global[ 159 ]++;
      if(nUnused == 0 && myfeof(zStream)) break;
   }
   clava_dcg_global[ 153 ]++;
   if(ferror(zStream)) goto errhandler_io;
   clava_dcg_global[ 160 ]++;
   ret = fclose(zStream);
   if(ret == (-1)) goto errhandler_io;
   if(verbosity >= 2) {
      clava_dcg_global[ 161 ]++;
      fprintf(stderr, "\n    ");
   }
   
   return ((Bool) 1);
   errhandler:
   clava_dcg_global[ 158 ]++;
   BZ2_bzReadClose(&bzerr_dummy, bzf);
   if(verbosity == 0) {
      clava_dcg_global[ 161 ]++;
      fprintf(stderr, "%s: %s: ", progName, inName);
   }
   switch (bzerr) {
      case (-9):
      clava_dcg_global[ 162 ]++;
      configError();
      break;
      case (-6):
      errhandler_io:
      clava_dcg_global[ 163 ]++;
      ioError();
      break;
      case (-4):
      clava_dcg_global[ 161 ]++;
      fprintf(stderr, "data integrity (CRC) error in data\n");
      
      return ((Bool) 0);
      case (-3):
      clava_dcg_global[ 164 ]++;
      outOfMemory();
      case (-7):
      clava_dcg_global[ 161 ]++;
      fprintf(stderr, "file ends unexpectedly\n");
      
      return ((Bool) 0);
      case (-5):
      if(zStream != stdin) {
         clava_dcg_global[ 160 ]++;
         fclose(zStream);
      }
      if(streamNo == 1) {
         clava_dcg_global[ 161 ]++;
         fprintf(stderr, "bad magic number (file not created by bzip2)\n");
         
         return ((Bool) 0);
      }
      else {
         if(noisy) {
            clava_dcg_global[ 161 ]++;
            fprintf(stderr, "trailing garbage after EOF ignored\n");
         }
         
         return ((Bool) 1);
      }
      default:
      clava_dcg_global[ 157 ]++;
      panic("test:unexpected error");
   }
   clava_dcg_global[ 157 ]++;
   panic("test:end");
   
   return ((Bool) 1);
   /*notreached*/
}

/*---------------------------------------------------*/
/*--- Error [non-] handling grunge                ---*/
/*---------------------------------------------------*/
/*---------------------------------------------*/
static void setExit(Int32 v) {
   if(v > exitValue) exitValue = v;
}

/*---------------------------------------------*/
static void cadvise() {
   if(noisy) {
      clava_dcg_global[ 165 ]++;
      fprintf(stderr, "\nIt is possible that the compressed file(s) have become corrupted.\nYou can use the -tvv option to test integrity of such files.\n\nYou can use the `bzip2recover' program to attempt to recover\ndata from undamaged sections of corrupted files.\n\n");
   }
}

/*---------------------------------------------*/
static void showFileNames() {
   if(noisy) {
      clava_dcg_global[ 166 ]++;
      fprintf(stderr, "\tInput file = %s, output file = %s\n", inName, outName);
   }
}

/*---------------------------------------------*/
static void cleanUpAndFail(Int32 ec) {
   IntNative retVal;
   struct stat statBuf;
   if(srcMode == 3 && opMode != 3 && deleteOutputOnInterrupt) {
      /*Check whether input file still exists.  Delete output file
      only if input exists to avoid loss of data.  Joerg Prante, 5
      January 2002.  (JRS 06-Jan-2002: other changes in 1.0.2 mean
      this is less likely to happen.  But to be ultra-paranoid, we
      do the check anyway.)*/
      clava_dcg_global[ 167 ]++;
      retVal = stat(inName, &statBuf);
      if(retVal == 0) {
         if(noisy) {
            clava_dcg_global[ 168 ]++;
            fprintf(stderr, "%s: Deleting output file %s, if it exists.\n", progName, outName);
         }
         if(outputHandleJustInCase != ((void *) 0)) {
            clava_dcg_global[ 169 ]++;
            fclose(outputHandleJustInCase);
         }
         clava_dcg_global[ 170 ]++;
         retVal = remove(outName);
         if(retVal != 0) {
            clava_dcg_global[ 168 ]++;
            fprintf(stderr, "%s: WARNING: deletion of output file (apparently) failed.\n", progName);
         }
      }
      else {
         clava_dcg_global[ 168 ]++;
         fprintf(stderr, "%s: WARNING: deletion of output file suppressed\n", progName);
         clava_dcg_global[ 168 ]++;
         fprintf(stderr, "%s:    since input file no longer exists.  Output file\n", progName);
         clava_dcg_global[ 168 ]++;
         fprintf(stderr, "%s:    `%s' may be incomplete.\n", progName, outName);
         clava_dcg_global[ 168 ]++;
         fprintf(stderr, "%s:    I suggest doing an integrity test (bzip2 -tv) of it.\n", progName);
      }
   }
   if(noisy && numFileNames > 0 && numFilesProcessed < numFileNames) {
      clava_dcg_global[ 168 ]++;
      fprintf(stderr, "%s: WARNING: some files have not been processed:\n%s:    %d specified on command line, %d not processed yet.\n\n", progName, progName, numFileNames, numFileNames - numFilesProcessed);
   }
   clava_dcg_global[ 171 ]++;
   setExit(ec);
   clava_dcg_global[ 172 ]++;
   exit(exitValue);
}

/*---------------------------------------------*/
static void panic(Char *s) {
   clava_dcg_global[ 173 ]++;
   fprintf(stderr, "\n%s: PANIC -- internal consistency error:\n\t%s\n\tThis is a BUG.  Please report it to me at:\n\tjseward@acm.org\n", progName, s);
   clava_dcg_global[ 174 ]++;
   showFileNames();
   clava_dcg_global[ 175 ]++;
   cleanUpAndFail(3);
}

/*---------------------------------------------*/
static void crcError() {
   clava_dcg_global[ 176 ]++;
   fprintf(stderr, "\n%s: Data integrity error when decompressing.\n", progName);
   clava_dcg_global[ 177 ]++;
   showFileNames();
   clava_dcg_global[ 178 ]++;
   cadvise();
   clava_dcg_global[ 179 ]++;
   cleanUpAndFail(2);
}

/*---------------------------------------------*/
static void compressedStreamEOF() {
   if(noisy) {
      clava_dcg_global[ 180 ]++;
      fprintf(stderr, "\n%s: Compressed file ends unexpectedly;\n\tperhaps it is corrupted?  *Possible* reason follows.\n", progName);
      clava_dcg_global[ 181 ]++;
      perror(progName);
      clava_dcg_global[ 182 ]++;
      showFileNames();
      clava_dcg_global[ 183 ]++;
      cadvise();
   }
   clava_dcg_global[ 184 ]++;
   cleanUpAndFail(2);
}

/*---------------------------------------------*/
static void ioError() {
   clava_dcg_global[ 185 ]++;
   fprintf(stderr, "\n%s: I/O or other error, bailing out.  Possible reason follows.\n", progName);
   clava_dcg_global[ 186 ]++;
   perror(progName);
   clava_dcg_global[ 187 ]++;
   showFileNames();
   clava_dcg_global[ 188 ]++;
   cleanUpAndFail(1);
}

/*---------------------------------------------*/
static void mySignalCatcher(IntNative n) {
   clava_dcg_global[ 189 ]++;
   fprintf(stderr, "\n%s: Control-C or similar caught, quitting.\n", progName);
   clava_dcg_global[ 190 ]++;
   cleanUpAndFail(1);
}

/*---------------------------------------------*/
static void mySIGSEGVorSIGBUScatcher(IntNative n) {
   if(opMode == 1) {
      clava_dcg_global[ 191 ]++;
      fprintf(stderr, "\n%s: Caught a SIGSEGV or SIGBUS whilst compressing.\n\n   Possible causes are (most likely first):\n   (1) This computer has unreliable memory or cache hardware\n       (a surprisingly common problem; try a different machine.)\n   (2) A bug in the compiler used to create this executable\n       (unlikely, if you didn't compile bzip2 yourself.)\n   (3) A real bug in bzip2 -- I hope this should never be the case.\n   The user's manual, Section 4.3, has more info on (1) and (2).\n   \n   If you suspect this is a bug in bzip2, or are unsure about (1)\n   or (2), feel free to report it to me at: jseward@acm.org.\n   Section 4.3 of the user's manual describes the info a useful\n   bug report should have.  If the manual is available on your\n   system, please try and read it before mailing me.  If you don't\n   have the manual or can't be bothered to read it, mail me anyway.\n\n", progName);
   }
   else {
      clava_dcg_global[ 191 ]++;
      fprintf(stderr, "\n%s: Caught a SIGSEGV or SIGBUS whilst decompressing.\n\n   Possible causes are (most likely first):\n   (1) The compressed data is corrupted, and bzip2's usual checks\n       failed to detect this.  Try bzip2 -tvv my_file.bz2.\n   (2) This computer has unreliable memory or cache hardware\n       (a surprisingly common problem; try a different machine.)\n   (3) A bug in the compiler used to create this executable\n       (unlikely, if you didn't compile bzip2 yourself.)\n   (4) A real bug in bzip2 -- I hope this should never be the case.\n   The user's manual, Section 4.3, has more info on (2) and (3).\n   \n   If you suspect this is a bug in bzip2, or are unsure about (2)\n   or (3), feel free to report it to me at: jseward@acm.org.\n   Section 4.3 of the user's manual describes the info a useful\n   bug report should have.  If the manual is available on your\n   system, please try and read it before mailing me.  If you don't\n   have the manual or can't be bothered to read it, mail me anyway.\n\n", progName);
   }
   clava_dcg_global[ 192 ]++;
   showFileNames();
   if(opMode == 1) {
      clava_dcg_global[ 193 ]++;
      cleanUpAndFail(3);
   }
   else {
      clava_dcg_global[ 194 ]++;
      cadvise();
      clava_dcg_global[ 193 ]++;
      cleanUpAndFail(2);
   }
}

/*---------------------------------------------*/
static void outOfMemory() {
   clava_dcg_global[ 195 ]++;
   fprintf(stderr, "\n%s: couldn't allocate enough memory\n", progName);
   clava_dcg_global[ 196 ]++;
   showFileNames();
   clava_dcg_global[ 197 ]++;
   cleanUpAndFail(1);
}

/*---------------------------------------------*/
static void configError() {
   clava_dcg_global[ 198 ]++;
   fprintf(stderr, "bzip2: I'm not configured correctly for this platform!\n\tI require Int32, Int16 and Char to have sizes\n\tof 4, 2 and 1 bytes to run properly, and they don't.\n\tProbably you can fix this by defining them correctly,\n\tand recompiling.  Bye!\n");
   clava_dcg_global[ 199 ]++;
   setExit(3);
   clava_dcg_global[ 200 ]++;
   exit(exitValue);
}

/*---------------------------------------------------*/
/*--- The main driver machinery                   ---*/
/*---------------------------------------------------*/
/*All rather crufty.  The main problem is that input files
are stat()d multiple times before use.  This should be
cleaned up.
*/
/*---------------------------------------------*/
static void pad(Char *s) {
   Int32 i;
   clava_dcg_global[ 201 ]++;
   if((Int32) strlen(s) >= longestFileName) 
   return;
   for(i = 1; i <= longestFileName - (Int32) strlen(s); i++) {
      clava_dcg_global[ 201 ]++;
      clava_dcg_global[ 202 ]++;
      fprintf(stderr, " ");
   }
}

/*---------------------------------------------*/
static void copyFileName(Char *to, Char *from) {
   clava_dcg_global[ 203 ]++;
   if(strlen(from) > 1034 - 10) {
      clava_dcg_global[ 204 ]++;
      fprintf(stderr, "bzip2: file name\n`%s'\nis suspiciously (more than %d chars) long.\nTry using a reasonable file name instead.  Sorry! :-)\n", from, 1034 - 10);
      clava_dcg_global[ 205 ]++;
      setExit(1);
      clava_dcg_global[ 206 ]++;
      exit(exitValue);
   }
   clava_dcg_global[ 207 ]++;
   strncpy(to, from, 1034 - 10);
   to[1034 - 10] = '\0';
}

/*---------------------------------------------*/
static Bool fileExists(Char *name) {
   clava_dcg_global[ 208 ]++;
   FILE *tmp = fopen(name, "rb");
   Bool exists = (tmp != ((void *) 0));
   if(tmp != ((void *) 0)) {
      clava_dcg_global[ 209 ]++;
      fclose(tmp);
   }
   
   return exists;
}

/*---------------------------------------------*/
/*Open an output file safely with O_EXCL and good permissions.
This avoids a race condition in versions < 1.0.2, in which
the file was first opened and then had its interim permissions
set safely.  We instead use open() to create the file with
the interim permissions required. (--- --- rw-).

For non-Unix platforms, if we are not worrying about
security issues, simple this simply behaves like fopen.
*/
FILE * fopen_output_safely(Char *name, char const *mode) {
   FILE *fp;
   IntNative fh;
   clava_dcg_global[ 210 ]++;
   fh = open(name, 01 | 0100 | 0200, 0200 | 0400);
   if(fh == -1) 
   return ((void *) 0);
   clava_dcg_global[ 211 ]++;
   fp = fdopen(fh, mode);
   if(fp == ((void *) 0)) {
      clava_dcg_global[ 212 ]++;
      close(fh);
   }
   
   return fp;
}

/*---------------------------------------------*/
/*--
if in doubt, return True
--*/
static Bool notAStandardFile(Char *name) {
   IntNative i;
   struct stat statBuf;
   clava_dcg_global[ 213 ]++;
   i = lstat(name, &statBuf);
   if(i != 0) 
   return ((Bool) 1);
   if(((((statBuf.st_mode)) & 0170000) == (0100000))) 
   return ((Bool) 0);
   
   return ((Bool) 1);
}

/*---------------------------------------------*/
/*--
rac 11/21/98 see if file has hard links to it
--*/
static Int32 countHardLinks(Char *name) {
   IntNative i;
   struct stat statBuf;
   clava_dcg_global[ 214 ]++;
   i = lstat(name, &statBuf);
   if(i != 0) 
   return 0;
   
   return (statBuf.st_nlink - 1);
}

/*---------------------------------------------*/

/*Copy modification date, access date, permissions and owner from the
source to destination file.  We have to copy this meta-info off
into fileMetaInfo before starting to compress / decompress it,
because doing it afterwards means we get the wrong access time.

To complicate matters, in compress() and decompress() below, the
sequence of tests preceding the call to saveInputFileMetaInfo()
involves calling fileExists(), which in turn establishes its result
by attempting to fopen() the file, and if successful, immediately
fclose()ing it again.  So we have to assume that the fopen() call
does not cause the access time field to be updated.

Reading of the man page for stat() (man 2 stat) on RedHat 7.2 seems
to imply that merely doing open() will not affect the access time.
Therefore we merely need to hope that the C library only does
open() as a result of fopen(), and not any kind of read()-ahead
cleverness.

It sounds pretty fragile to me.  Whether this carries across
robustly to arbitrary Unix-like platforms (or even works robustly
on this one, RedHat 7.2) is unknown to me.  Nevertheless ...
*/

static struct stat fileMetaInfo;
static void saveInputFileMetaInfo(Char *srcName) {
   IntNative retVal;
   /*Note use of stat here, not lstat.*/
   clava_dcg_global[ 215 ]++;
   retVal = stat(srcName, &fileMetaInfo);
   {
      if((retVal) != 0) {
         clava_dcg_global[ 216 ]++;
         ioError();
      }
   }
   ;
}

static void applySavedMetaInfoToOutputFile(Char *dstName) {
   IntNative retVal;
   struct utimbuf uTimBuf;
   uTimBuf.actime = fileMetaInfo.st_atim.tv_sec;
   uTimBuf.modtime = fileMetaInfo.st_mtim.tv_sec;
   clava_dcg_global[ 217 ]++;
   retVal = chmod(dstName, fileMetaInfo.st_mode);
   {
      if((retVal) != 0) {
         clava_dcg_global[ 218 ]++;
         ioError();
      }
   }
   ;
   clava_dcg_global[ 219 ]++;
   retVal = utime(dstName, &uTimBuf);
   {
      if((retVal) != 0) {
         clava_dcg_global[ 218 ]++;
         ioError();
      }
   }
   ;
   clava_dcg_global[ 220 ]++;
   retVal = chown(dstName, fileMetaInfo.st_uid, fileMetaInfo.st_gid);
   /*chown() will in many cases return with EPERM, which can
   be safely ignored.
   */
}

/*---------------------------------------------*/
static Bool containsDubiousChars(Char *name) {
   /*On unix, files can contain any characters and the file expansion
   * is performed by the shell.
   */
   
   return ((Bool) 0);
   /*! BZ_UNIX*/
   /*On non-unix (Win* platforms), wildcard characters are not allowed in
   * filenames.
   */
   /*BZ_UNIX*/
}

/*---------------------------------------------*/

Char *zSuffix[4] = {".bz2", ".bz", ".tbz2", ".tbz"};
Char *unzSuffix[4] = {"", "", ".tar", ".tar"};
static Bool hasSuffix(Char *s, Char *suffix) {
   clava_dcg_global[ 221 ]++;
   Int32 ns = strlen(s);
   clava_dcg_global[ 221 ]++;
   Int32 nx = strlen(suffix);
   if(ns < nx) 
   return ((Bool) 0);
   clava_dcg_global[ 222 ]++;
   if(strcmp(s + ns - nx, suffix) == 0) 
   return ((Bool) 1);
   
   return ((Bool) 0);
}

static Bool mapSuffix(Char *name, Char *oldSuffix, Char *newSuffix) {
   clava_dcg_global[ 223 ]++;
   if(!hasSuffix(name, oldSuffix)) 
   return ((Bool) 0);
   clava_dcg_global[ 224 ]++;
   clava_dcg_global[ 224 ]++;
   name[strlen(name) - strlen(oldSuffix)] = 0;
   clava_dcg_global[ 225 ]++;
   strcat(name, newSuffix);
   
   return ((Bool) 1);
}

/*---------------------------------------------*/
static void compress(Char *name) {
   FILE *inStr;
   FILE *outStr;
   Int32 n, i;
   struct stat statBuf;
   deleteOutputOnInterrupt = ((Bool) 0);
   if(name == ((void *) 0) && srcMode != 1) {
      clava_dcg_global[ 226 ]++;
      panic("compress: bad modes\n");
   }
   switch (srcMode) {
      case 1:
      clava_dcg_global[ 227 ]++;
      copyFileName(inName, "(stdin)");
      clava_dcg_global[ 227 ]++;
      copyFileName(outName, "(stdout)");
      break;
      case 3:
      clava_dcg_global[ 227 ]++;
      copyFileName(inName, name);
      clava_dcg_global[ 227 ]++;
      copyFileName(outName, name);
      clava_dcg_global[ 228 ]++;
      strcat(outName, ".bz2");
      break;
      case 2:
      clava_dcg_global[ 227 ]++;
      copyFileName(inName, name);
      clava_dcg_global[ 227 ]++;
      copyFileName(outName, "(stdout)");
      break;
   }
   clava_dcg_global[ 229 ]++;
   if(srcMode != 1 && containsDubiousChars(inName)) {
      if(noisy) {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: There are no files matching `%s'.\n", progName, inName);
      }
      clava_dcg_global[ 231 ]++;
      setExit(1);
      
      return;
   }
   clava_dcg_global[ 232 ]++;
   if(srcMode != 1 && !fileExists(inName)) {
      clava_dcg_global[ 230 ]++;
      clava_dcg_global[ 233 ]++;
      clava_dcg_global[ 234 ]++;
      fprintf(stderr, "%s: Can't open input file %s: %s.\n", progName, inName, strerror((*__errno_location())));
      clava_dcg_global[ 231 ]++;
      setExit(1);
      
      return;
   }
   for(i = 0; i < 4; i++) {
      clava_dcg_global[ 235 ]++;
      if(hasSuffix(inName, zSuffix[i])) {
         if(noisy) {
            clava_dcg_global[ 230 ]++;
            fprintf(stderr, "%s: Input file %s already has %s suffix.\n", progName, inName, zSuffix[i]);
         }
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
   }
   if(srcMode == 3 || srcMode == 2) {
      clava_dcg_global[ 236 ]++;
      stat(inName, &statBuf);
      if(((((statBuf.st_mode)) & 0170000) == (0040000))) {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: Input file %s is a directory.\n", progName, inName);
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
   }
   clava_dcg_global[ 237 ]++;
   if(srcMode == 3 && !forceOverwrite && notAStandardFile(inName)) {
      if(noisy) {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: Input file %s is not a normal file.\n", progName, inName);
      }
      clava_dcg_global[ 231 ]++;
      setExit(1);
      
      return;
   }
   clava_dcg_global[ 232 ]++;
   if(srcMode == 3 && fileExists(outName)) {
      if(forceOverwrite) {
         clava_dcg_global[ 238 ]++;
         remove(outName);
      }
      else {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: Output file %s already exists.\n", progName, outName);
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
   }
   clava_dcg_global[ 239 ]++;
   if(srcMode == 3 && !forceOverwrite && (n = countHardLinks(inName)) > 0) {
      clava_dcg_global[ 230 ]++;
      fprintf(stderr, "%s: Input file %s has %d other link%s.\n", progName, inName, n, n > 1 ? "s" : "");
      clava_dcg_global[ 231 ]++;
      setExit(1);
      
      return;
   }
   if(srcMode == 3) {
      /*Save the file's meta-info before we open it.  Doing it later
      means we mess up the access times.*/
      clava_dcg_global[ 240 ]++;
      saveInputFileMetaInfo(inName);
   }
   switch (srcMode) {
      case 1:
      inStr = stdin;
      outStr = stdout;
      clava_dcg_global[ 241 ]++;
      clava_dcg_global[ 242 ]++;
      if(isatty(fileno(stdout))) {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: I won't write compressed data to a terminal.\n", progName);
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: For help, type: `%s --help'.\n", progName, progName);
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      case 2:
      clava_dcg_global[ 243 ]++;
      inStr = fopen(inName, "rb");
      outStr = stdout;
      clava_dcg_global[ 241 ]++;
      clava_dcg_global[ 242 ]++;
      if(isatty(fileno(stdout))) {
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: I won't write compressed data to a terminal.\n", progName);
         clava_dcg_global[ 230 ]++;
         fprintf(stderr, "%s: For help, type: `%s --help'.\n", progName, progName);
         if(inStr != ((void *) 0)) {
            clava_dcg_global[ 244 ]++;
            fclose(inStr);
         }
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
      ;
      if(inStr == ((void *) 0)) {
         clava_dcg_global[ 230 ]++;
         clava_dcg_global[ 233 ]++;
         clava_dcg_global[ 234 ]++;
         fprintf(stderr, "%s: Can't open input file %s: %s.\n", progName, inName, strerror((*__errno_location())));
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      case 3:
      clava_dcg_global[ 243 ]++;
      inStr = fopen(inName, "rb");
      clava_dcg_global[ 245 ]++;
      outStr = fopen_output_safely(outName, "wb");
      if(outStr == ((void *) 0)) {
         clava_dcg_global[ 230 ]++;
         clava_dcg_global[ 233 ]++;
         clava_dcg_global[ 234 ]++;
         fprintf(stderr, "%s: Can't create output file %s: %s.\n", progName, outName, strerror((*__errno_location())));
         if(inStr != ((void *) 0)) {
            clava_dcg_global[ 244 ]++;
            fclose(inStr);
         }
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
      if(inStr == ((void *) 0)) {
         clava_dcg_global[ 230 ]++;
         clava_dcg_global[ 233 ]++;
         clava_dcg_global[ 234 ]++;
         fprintf(stderr, "%s: Can't open input file %s: %s.\n", progName, inName, strerror((*__errno_location())));
         if(outStr != ((void *) 0)) {
            clava_dcg_global[ 244 ]++;
            fclose(outStr);
         }
         clava_dcg_global[ 231 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      default:
      clava_dcg_global[ 226 ]++;
      panic("compress: bad srcMode");
      break;
   }
   if(verbosity >= 1) {
      clava_dcg_global[ 230 ]++;
      fprintf(stderr, "  %s: ", inName);
      clava_dcg_global[ 246 ]++;
      pad(inName);
      clava_dcg_global[ 247 ]++;
      fflush(stderr);
   }
   /*--- Now the input and output handles are sane.  Do the Biz. ---*/
   outputHandleJustInCase = outStr;
   deleteOutputOnInterrupt = ((Bool) 1);
   clava_dcg_global[ 248 ]++;
   compressStream(inStr, outStr);
   outputHandleJustInCase = ((void *) 0);
   /*--- If there was an I/O error, we won't get here. ---*/
   if(srcMode == 3) {
      clava_dcg_global[ 249 ]++;
      applySavedMetaInfoToOutputFile(outName);
      deleteOutputOnInterrupt = ((Bool) 0);
      if(!keepInputFiles) {
         clava_dcg_global[ 238 ]++;
         IntNative retVal = remove(inName);
         {
            if((retVal) != 0) {
               clava_dcg_global[ 250 ]++;
               ioError();
            }
         }
         ;
      }
   }
   deleteOutputOnInterrupt = ((Bool) 0);
}

/*---------------------------------------------*/
static void uncompress(Char *name) {
   FILE *inStr;
   FILE *outStr;
   Int32 n, i;
   Bool magicNumberOK;
   Bool cantGuess;
   struct stat statBuf;
   deleteOutputOnInterrupt = ((Bool) 0);
   if(name == ((void *) 0) && srcMode != 1) {
      clava_dcg_global[ 251 ]++;
      panic("uncompress: bad modes\n");
   }
   cantGuess = ((Bool) 0);
   switch (srcMode) {
      case 1:
      clava_dcg_global[ 252 ]++;
      copyFileName(inName, "(stdin)");
      clava_dcg_global[ 252 ]++;
      copyFileName(outName, "(stdout)");
      break;
      case 3:
      clava_dcg_global[ 252 ]++;
      copyFileName(inName, name);
      clava_dcg_global[ 252 ]++;
      copyFileName(outName, name);
      for(i = 0; i < 4; i++) {
         clava_dcg_global[ 253 ]++;
         if(mapSuffix(outName, zSuffix[i], unzSuffix[i])) goto zzz;
      }
      cantGuess = ((Bool) 1);
      clava_dcg_global[ 254 ]++;
      strcat(outName, ".out");
      break;
      case 2:
      clava_dcg_global[ 252 ]++;
      copyFileName(inName, name);
      clava_dcg_global[ 252 ]++;
      copyFileName(outName, "(stdout)");
      break;
   }
   zzz:
   clava_dcg_global[ 255 ]++;
   if(srcMode != 1 && containsDubiousChars(inName)) {
      if(noisy) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: There are no files matching `%s'.\n", progName, inName);
      }
      clava_dcg_global[ 257 ]++;
      setExit(1);
      
      return;
   }
   clava_dcg_global[ 258 ]++;
   if(srcMode != 1 && !fileExists(inName)) {
      clava_dcg_global[ 256 ]++;
      clava_dcg_global[ 259 ]++;
      clava_dcg_global[ 260 ]++;
      fprintf(stderr, "%s: Can't open input file %s: %s.\n", progName, inName, strerror((*__errno_location())));
      clava_dcg_global[ 257 ]++;
      setExit(1);
      
      return;
   }
   if(srcMode == 3 || srcMode == 2) {
      clava_dcg_global[ 261 ]++;
      stat(inName, &statBuf);
      if(((((statBuf.st_mode)) & 0170000) == (0040000))) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: Input file %s is a directory.\n", progName, inName);
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
   }
   clava_dcg_global[ 262 ]++;
   if(srcMode == 3 && !forceOverwrite && notAStandardFile(inName)) {
      if(noisy) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: Input file %s is not a normal file.\n", progName, inName);
      }
      clava_dcg_global[ 257 ]++;
      setExit(1);
      
      return;
   }
   if(cantGuess) {
      /*srcMode == SM_F2F implied &&*/
      if(noisy) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: Can't guess original name for %s -- using %s\n", progName, inName, outName);
      }
      /*just a warning, no return*/
   }
   clava_dcg_global[ 258 ]++;
   if(srcMode == 3 && fileExists(outName)) {
      if(forceOverwrite) {
         clava_dcg_global[ 263 ]++;
         remove(outName);
      }
      else {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: Output file %s already exists.\n", progName, outName);
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
   }
   clava_dcg_global[ 264 ]++;
   if(srcMode == 3 && !forceOverwrite && (n = countHardLinks(inName)) > 0) {
      clava_dcg_global[ 256 ]++;
      fprintf(stderr, "%s: Input file %s has %d other link%s.\n", progName, inName, n, n > 1 ? "s" : "");
      clava_dcg_global[ 257 ]++;
      setExit(1);
      
      return;
   }
   if(srcMode == 3) {
      /*Save the file's meta-info before we open it.  Doing it later
      means we mess up the access times.*/
      clava_dcg_global[ 265 ]++;
      saveInputFileMetaInfo(inName);
   }
   switch (srcMode) {
      case 1:
      inStr = stdin;
      outStr = stdout;
      clava_dcg_global[ 266 ]++;
      clava_dcg_global[ 267 ]++;
      if(isatty(fileno(stdin))) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: I won't read compressed data from a terminal.\n", progName);
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: For help, type: `%s --help'.\n", progName, progName);
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      case 2:
      clava_dcg_global[ 268 ]++;
      inStr = fopen(inName, "rb");
      outStr = stdout;
      if(inStr == ((void *) 0)) {
         clava_dcg_global[ 256 ]++;
         clava_dcg_global[ 259 ]++;
         clava_dcg_global[ 260 ]++;
         fprintf(stderr, "%s: Can't open input file %s:%s.\n", progName, inName, strerror((*__errno_location())));
         if(inStr != ((void *) 0)) {
            clava_dcg_global[ 269 ]++;
            fclose(inStr);
         }
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      case 3:
      clava_dcg_global[ 268 ]++;
      inStr = fopen(inName, "rb");
      clava_dcg_global[ 270 ]++;
      outStr = fopen_output_safely(outName, "wb");
      if(outStr == ((void *) 0)) {
         clava_dcg_global[ 256 ]++;
         clava_dcg_global[ 259 ]++;
         clava_dcg_global[ 260 ]++;
         fprintf(stderr, "%s: Can't create output file %s: %s.\n", progName, outName, strerror((*__errno_location())));
         if(inStr != ((void *) 0)) {
            clava_dcg_global[ 269 ]++;
            fclose(inStr);
         }
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
      if(inStr == ((void *) 0)) {
         clava_dcg_global[ 256 ]++;
         clava_dcg_global[ 259 ]++;
         clava_dcg_global[ 260 ]++;
         fprintf(stderr, "%s: Can't open input file %s: %s.\n", progName, inName, strerror((*__errno_location())));
         if(outStr != ((void *) 0)) {
            clava_dcg_global[ 269 ]++;
            fclose(outStr);
         }
         clava_dcg_global[ 257 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      default:
      clava_dcg_global[ 251 ]++;
      panic("uncompress: bad srcMode");
      break;
   }
   if(verbosity >= 1) {
      clava_dcg_global[ 256 ]++;
      fprintf(stderr, "  %s: ", inName);
      clava_dcg_global[ 271 ]++;
      pad(inName);
      clava_dcg_global[ 272 ]++;
      fflush(stderr);
   }
   /*--- Now the input and output handles are sane.  Do the Biz. ---*/
   outputHandleJustInCase = outStr;
   deleteOutputOnInterrupt = ((Bool) 1);
   clava_dcg_global[ 273 ]++;
   magicNumberOK = uncompressStream(inStr, outStr);
   outputHandleJustInCase = ((void *) 0);
   /*--- If there was an I/O error, we won't get here. ---*/
   if(magicNumberOK) {
      if(srcMode == 3) {
         clava_dcg_global[ 274 ]++;
         applySavedMetaInfoToOutputFile(outName);
         deleteOutputOnInterrupt = ((Bool) 0);
         if(!keepInputFiles) {
            clava_dcg_global[ 263 ]++;
            IntNative retVal = remove(inName);
            {
               if((retVal) != 0) {
                  clava_dcg_global[ 275 ]++;
                  ioError();
               }
            }
            ;
         }
      }
   }
   else {
      unzFailsExist = ((Bool) 1);
      deleteOutputOnInterrupt = ((Bool) 0);
      if(srcMode == 3) {
         clava_dcg_global[ 263 ]++;
         IntNative retVal = remove(outName);
         {
            if((retVal) != 0) {
               clava_dcg_global[ 275 ]++;
               ioError();
            }
         }
         ;
      }
   }
   deleteOutputOnInterrupt = ((Bool) 0);
   if(magicNumberOK) {
      if(verbosity >= 1) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "done\n");
      }
   }
   else {
      clava_dcg_global[ 257 ]++;
      setExit(2);
      if(verbosity >= 1) {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "not a bzip2 file.\n");
      }
      else {
         clava_dcg_global[ 256 ]++;
         fprintf(stderr, "%s: %s is not a bzip2 file.\n", progName, inName);
      }
   }
}

/*---------------------------------------------*/
static void testf(Char *name) {
   FILE *inStr;
   Bool allOK;
   struct stat statBuf;
   deleteOutputOnInterrupt = ((Bool) 0);
   if(name == ((void *) 0) && srcMode != 1) {
      clava_dcg_global[ 276 ]++;
      panic("testf: bad modes\n");
   }
   clava_dcg_global[ 277 ]++;
   copyFileName(outName, "(none)");
   switch (srcMode) {
      case 1:
      clava_dcg_global[ 277 ]++;
      copyFileName(inName, "(stdin)");
      break;
      case 3:
      clava_dcg_global[ 277 ]++;
      copyFileName(inName, name);
      break;
      case 2:
      clava_dcg_global[ 277 ]++;
      copyFileName(inName, name);
      break;
   }
   clava_dcg_global[ 278 ]++;
   if(srcMode != 1 && containsDubiousChars(inName)) {
      if(noisy) {
         clava_dcg_global[ 279 ]++;
         fprintf(stderr, "%s: There are no files matching `%s'.\n", progName, inName);
      }
      clava_dcg_global[ 280 ]++;
      setExit(1);
      
      return;
   }
   clava_dcg_global[ 281 ]++;
   if(srcMode != 1 && !fileExists(inName)) {
      clava_dcg_global[ 279 ]++;
      clava_dcg_global[ 282 ]++;
      clava_dcg_global[ 283 ]++;
      fprintf(stderr, "%s: Can't open input %s: %s.\n", progName, inName, strerror((*__errno_location())));
      clava_dcg_global[ 280 ]++;
      setExit(1);
      
      return;
   }
   if(srcMode != 1) {
      clava_dcg_global[ 284 ]++;
      stat(inName, &statBuf);
      if(((((statBuf.st_mode)) & 0170000) == (0040000))) {
         clava_dcg_global[ 279 ]++;
         fprintf(stderr, "%s: Input file %s is a directory.\n", progName, inName);
         clava_dcg_global[ 280 ]++;
         setExit(1);
         
         return;
      }
   }
   switch (srcMode) {
      case 1:
      clava_dcg_global[ 285 ]++;
      clava_dcg_global[ 286 ]++;
      if(isatty(fileno(stdin))) {
         clava_dcg_global[ 279 ]++;
         fprintf(stderr, "%s: I won't read compressed data from a terminal.\n", progName);
         clava_dcg_global[ 279 ]++;
         fprintf(stderr, "%s: For help, type: `%s --help'.\n", progName, progName);
         clava_dcg_global[ 280 ]++;
         setExit(1);
         
         return;
      }
      ;
      inStr = stdin;
      break;
      case 2:
      case 3:
      clava_dcg_global[ 287 ]++;
      inStr = fopen(inName, "rb");
      if(inStr == ((void *) 0)) {
         clava_dcg_global[ 279 ]++;
         clava_dcg_global[ 282 ]++;
         clava_dcg_global[ 283 ]++;
         fprintf(stderr, "%s: Can't open input file %s:%s.\n", progName, inName, strerror((*__errno_location())));
         clava_dcg_global[ 280 ]++;
         setExit(1);
         
         return;
      }
      ;
      break;
      default:
      clava_dcg_global[ 276 ]++;
      panic("testf: bad srcMode");
      break;
   }
   if(verbosity >= 1) {
      clava_dcg_global[ 279 ]++;
      fprintf(stderr, "  %s: ", inName);
      clava_dcg_global[ 288 ]++;
      pad(inName);
      clava_dcg_global[ 289 ]++;
      fflush(stderr);
   }
   /*--- Now the input handle is sane.  Do the Biz. ---*/
   outputHandleJustInCase = ((void *) 0);
   clava_dcg_global[ 290 ]++;
   allOK = testStream(inStr);
   if(allOK && verbosity >= 1) {
      clava_dcg_global[ 279 ]++;
      fprintf(stderr, "ok\n");
   }
   if(!allOK) testFailsExist = ((Bool) 1);
}

/*---------------------------------------------*/
static void license() {
   clava_dcg_global[ 291 ]++;
   clava_dcg_global[ 292 ]++;
   fprintf(stderr, "bzip2, a block-sorting file compressor.  Version %s.\n   \n   Copyright (C) 1996-2002 by Julian Seward.\n   \n   This program is free software; you can redistribute it and/or modify\n   it under the terms set out in the LICENSE file, which is included\n   in the bzip2-1.0 source distribution.\n   \n   This program is distributed in the hope that it will be useful,\n   but WITHOUT ANY WARRANTY; without even the implied warranty of\n   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n   LICENSE file for more details.\n   \n", BZ2_bzlibVersion());
}

/*---------------------------------------------*/
static void usage(Char *fullProgName) {
   clava_dcg_global[ 293 ]++;
   clava_dcg_global[ 294 ]++;
   fprintf(stderr, "bzip2, a block-sorting file compressor.  Version %s.\n\n   usage: %s [flags and input files in any order]\n\n   -h --help           print this message\n   -d --decompress     force decompression\n   -z --compress       force compression\n   -k --keep           keep (don't delete) input files\n   -f --force          overwrite existing output files\n   -t --test           test compressed file integrity\n   -c --stdout         output to standard out\n   -q --quiet          suppress noncritical error messages\n   -v --verbose        be verbose (a 2nd -v gives more)\n   -L --license        display software version & license\n   -V --version        display software version & license\n   -s --small          use less memory (at most 2500k)\n   -1 .. -9            set block size to 100k .. 900k\n   --fast              alias for -1\n   --best              alias for -9\n\n   If invoked as `bzip2', default action is to compress.\n              as `bunzip2',  default action is to decompress.\n              as `bzcat', default action is to decompress to stdout.\n\n   If no file names are given, bzip2 compresses or decompresses\n   from standard input to standard output.  You can combine\n   short flags, so `-v -4' means the same as -v4 or -4v, &c.\n\n", BZ2_bzlibVersion(), fullProgName);
}

/*---------------------------------------------*/
static void redundant(Char *flag) {
   clava_dcg_global[ 295 ]++;
   fprintf(stderr, "%s: %s is redundant in versions 0.9.5 and above\n", progName, flag);
}

/*---------------------------------------------*/
/*--
All the garbage from here to main() is purely to
implement a linked list of command-line arguments,
into which main() copies argv[1 .. argc-1].

The purpose of this exercise is to facilitate
the expansion of wildcard characters * and ? in
filenames for OSs which don't know how to do it
themselves, like MSDOS, Windows 95 and NT.

The actual Dirty Work is done by the platform-
specific macro APPEND_FILESPEC.
--*/

struct zzzz {
   Char *name;
   struct zzzz *link;
};

typedef struct zzzz Cell;
/*---------------------------------------------*/
static void * myMalloc(Int32 n) {
   void *p;
   clava_dcg_global[ 296 ]++;
   p = malloc((size_t) n);
   if(p == ((void *) 0)) {
      clava_dcg_global[ 297 ]++;
      outOfMemory();
   }
   
   return p;
}

/*---------------------------------------------*/
static Cell * mkCell() {
   Cell *c;
   clava_dcg_global[ 298 ]++;
   c = (Cell *) myMalloc(sizeof(Cell));
   c->name = ((void *) 0);
   c->link = ((void *) 0);
   
   return c;
}

/*---------------------------------------------*/
static Cell * snocString(Cell *root, Char *name) {
   if(root == ((void *) 0)) {
      clava_dcg_global[ 299 ]++;
      Cell *tmp = mkCell();
      clava_dcg_global[ 300 ]++;
      clava_dcg_global[ 301 ]++;
      tmp->name = (Char *) myMalloc(5 + strlen(name));
      clava_dcg_global[ 302 ]++;
      strcpy(tmp->name, name);
      
      return tmp;
   }
   else {
      Cell *tmp = root;
      while(tmp->link != ((void *) 0)) tmp = tmp->link;
      clava_dcg_global[ 303 ]++;
      tmp->link = snocString(tmp->link, name);
      
      return root;
   }
}

/*---------------------------------------------*/
static void addFlagsFromEnvVar(Cell **argList, Char *varName) {
   Int32 i, j, k;
   Char *envbase, *p;
   clava_dcg_global[ 304 ]++;
   envbase = getenv(varName);
   if(envbase != ((void *) 0)) {
      p = envbase;
      i = 0;
      while(((Bool) 1)) {
         if(p[i] == 0) break;
         p += i;
         i = 0;
         while(((*__ctype_b_loc())[(int) (((Int32) (p[0])))] & (unsigned short) _ISspace)) {
            clava_dcg_global[ 305 ]++;
            p++;
         }
         while(p[i] != 0 && !((*__ctype_b_loc())[(int) (((Int32) (p[i])))] & (unsigned short) _ISspace)) {
            clava_dcg_global[ 305 ]++;
            i++;
         }
         if(i > 0) {
            k = i;
            if(k > 1034 - 10) k = 1034 - 10;
            for(j = 0; j < k; j++) tmpName[j] = p[j];
            tmpName[k] = 0;
            clava_dcg_global[ 306 ]++;
            *argList = snocString((*argList), (tmpName));
         }
      }
   }
}

/*---------------------------------------------*/
void clava_call_graph() {
   FILE *log_file_365 = fopen("/home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/bzip2.dot", "w+");
   if (log_file_365 == NULL)
   {
       printf("Error opening file /home/specs/jbispo/repos/clava-examples/2019-C_Stress_Test/output/bzip2.dot\n");
       exit(1);
   } 
   fprintf(log_file_365, "digraph dynamic_call_graph {\n\n");
   if(clava_dcg_global[0] != 0) {
      fprintf(log_file_365, "	fallbackQSort3 -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[0]);
   }
   if(clava_dcg_global[1] != 0) {
      fprintf(log_file_365, "	fallbackQSort3 -> fallbackSimpleSort [label=\"%d\"];\n", clava_dcg_global[1]);
   }
   if(clava_dcg_global[2] != 0) {
      fprintf(log_file_365, "	fallbackSort -> fprintf [label=\"%d\"];\n", clava_dcg_global[2]);
   }
   if(clava_dcg_global[3] != 0) {
      fprintf(log_file_365, "	fallbackSort -> fallbackQSort3 [label=\"%d\"];\n", clava_dcg_global[3]);
   }
   if(clava_dcg_global[4] != 0) {
      fprintf(log_file_365, "	fallbackSort -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[4]);
   }
   if(clava_dcg_global[5] != 0) {
      fprintf(log_file_365, "	mainSimpleSort -> mainGtU [label=\"%d\"];\n", clava_dcg_global[5]);
   }
   if(clava_dcg_global[6] != 0) {
      fprintf(log_file_365, "	mainQSort3 -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[6]);
   }
   if(clava_dcg_global[7] != 0) {
      fprintf(log_file_365, "	mainQSort3 -> mainSimpleSort [label=\"%d\"];\n", clava_dcg_global[7]);
   }
   if(clava_dcg_global[8] != 0) {
      fprintf(log_file_365, "	mainQSort3 -> mmed3 [label=\"%d\"];\n", clava_dcg_global[8]);
   }
   if(clava_dcg_global[9] != 0) {
      fprintf(log_file_365, "	mainSort -> fprintf [label=\"%d\"];\n", clava_dcg_global[9]);
   }
   if(clava_dcg_global[10] != 0) {
      fprintf(log_file_365, "	mainSort -> mainQSort3 [label=\"%d\"];\n", clava_dcg_global[10]);
   }
   if(clava_dcg_global[11] != 0) {
      fprintf(log_file_365, "	mainSort -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[11]);
   }
   if(clava_dcg_global[12] != 0) {
      fprintf(log_file_365, "	BZ2_blockSort -> fallbackSort [label=\"%d\"];\n", clava_dcg_global[12]);
   }
   if(clava_dcg_global[13] != 0) {
      fprintf(log_file_365, "	BZ2_blockSort -> mainSort [label=\"%d\"];\n", clava_dcg_global[13]);
   }
   if(clava_dcg_global[14] != 0) {
      fprintf(log_file_365, "	BZ2_blockSort -> fprintf [label=\"%d\"];\n", clava_dcg_global[14]);
   }
   if(clava_dcg_global[15] != 0) {
      fprintf(log_file_365, "	BZ2_blockSort -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[15]);
   }
   if(clava_dcg_global[16] != 0) {
      fprintf(log_file_365, "	BZ2_hbMakeCodeLengths -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[16]);
   }
   if(clava_dcg_global[17] != 0) {
      fprintf(log_file_365, "	bsPutUInt32 -> bsW [label=\"%d\"];\n", clava_dcg_global[17]);
   }
   if(clava_dcg_global[18] != 0) {
      fprintf(log_file_365, "	bsPutUChar -> bsW [label=\"%d\"];\n", clava_dcg_global[18]);
   }
   if(clava_dcg_global[19] != 0) {
      fprintf(log_file_365, "	generateMTFValues -> makeMaps_e [label=\"%d\"];\n", clava_dcg_global[19]);
   }
   if(clava_dcg_global[20] != 0) {
      fprintf(log_file_365, "	sendMTFValues -> fprintf [label=\"%d\"];\n", clava_dcg_global[20]);
   }
   if(clava_dcg_global[21] != 0) {
      fprintf(log_file_365, "	sendMTFValues -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[21]);
   }
   if(clava_dcg_global[22] != 0) {
      fprintf(log_file_365, "	sendMTFValues -> BZ2_hbMakeCodeLengths [label=\"%d\"];\n", clava_dcg_global[22]);
   }
   if(clava_dcg_global[23] != 0) {
      fprintf(log_file_365, "	sendMTFValues -> BZ2_hbAssignCodes [label=\"%d\"];\n", clava_dcg_global[23]);
   }
   if(clava_dcg_global[24] != 0) {
      fprintf(log_file_365, "	sendMTFValues -> bsW [label=\"%d\"];\n", clava_dcg_global[24]);
   }
   if(clava_dcg_global[25] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> fprintf [label=\"%d\"];\n", clava_dcg_global[25]);
   }
   if(clava_dcg_global[26] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> BZ2_blockSort [label=\"%d\"];\n", clava_dcg_global[26]);
   }
   if(clava_dcg_global[27] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> BZ2_bsInitWrite [label=\"%d\"];\n", clava_dcg_global[27]);
   }
   if(clava_dcg_global[28] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> bsPutUChar [label=\"%d\"];\n", clava_dcg_global[28]);
   }
   if(clava_dcg_global[29] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> bsPutUInt32 [label=\"%d\"];\n", clava_dcg_global[29]);
   }
   if(clava_dcg_global[30] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> bsW [label=\"%d\"];\n", clava_dcg_global[30]);
   }
   if(clava_dcg_global[31] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> generateMTFValues [label=\"%d\"];\n", clava_dcg_global[31]);
   }
   if(clava_dcg_global[32] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> sendMTFValues [label=\"%d\"];\n", clava_dcg_global[32]);
   }
   if(clava_dcg_global[33] != 0) {
      fprintf(log_file_365, "	BZ2_compressBlock -> bsFinishWrite [label=\"%d\"];\n", clava_dcg_global[33]);
   }
   if(clava_dcg_global[34] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> strm [label=\"%d\"];\n", clava_dcg_global[34]);
   }
   if(clava_dcg_global[35] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> fprintf [label=\"%d\"];\n", clava_dcg_global[35]);
   }
   if(clava_dcg_global[36] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> makeMaps_d [label=\"%d\"];\n", clava_dcg_global[36]);
   }
   if(clava_dcg_global[37] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> BZ2_hbCreateDecodeTables [label=\"%d\"];\n", clava_dcg_global[37]);
   }
   if(clava_dcg_global[38] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> BZ2_indexIntoF [label=\"%d\"];\n", clava_dcg_global[38]);
   }
   if(clava_dcg_global[39] != 0) {
      fprintf(log_file_365, "	BZ2_decompress -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[39]);
   }
   if(clava_dcg_global[40] != 0) {
      fprintf(log_file_365, "	BZ2_bz__AssertH__fail -> fprintf [label=\"%d\"];\n", clava_dcg_global[40]);
   }
   if(clava_dcg_global[41] != 0) {
      fprintf(log_file_365, "	BZ2_bz__AssertH__fail -> BZ2_bzlibVersion [label=\"%d\"];\n", clava_dcg_global[41]);
   }
   if(clava_dcg_global[42] != 0) {
      fprintf(log_file_365, "	BZ2_bz__AssertH__fail -> exit [label=\"%d\"];\n", clava_dcg_global[42]);
   }
   if(clava_dcg_global[43] != 0) {
      fprintf(log_file_365, "	default_bzalloc -> malloc [label=\"%d\"];\n", clava_dcg_global[43]);
   }
   if(clava_dcg_global[44] != 0) {
      fprintf(log_file_365, "	default_bzfree -> free [label=\"%d\"];\n", clava_dcg_global[44]);
   }
   if(clava_dcg_global[45] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompressInit -> bz_config_ok [label=\"%d\"];\n", clava_dcg_global[45]);
   }
   if(clava_dcg_global[46] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompressInit -> strm [label=\"%d\"];\n", clava_dcg_global[46]);
   }
   if(clava_dcg_global[47] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompressInit -> init_RL [label=\"%d\"];\n", clava_dcg_global[47]);
   }
   if(clava_dcg_global[48] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompressInit -> prepare_new_block [label=\"%d\"];\n", clava_dcg_global[48]);
   }
   if(clava_dcg_global[49] != 0) {
      fprintf(log_file_365, "	flush_RL -> add_pair_to_block [label=\"%d\"];\n", clava_dcg_global[49]);
   }
   if(clava_dcg_global[50] != 0) {
      fprintf(log_file_365, "	flush_RL -> init_RL [label=\"%d\"];\n", clava_dcg_global[50]);
   }
   if(clava_dcg_global[51] != 0) {
      fprintf(log_file_365, "	copy_input_until_stop -> add_pair_to_block [label=\"%d\"];\n", clava_dcg_global[51]);
   }
   if(clava_dcg_global[52] != 0) {
      fprintf(log_file_365, "	handle_compress -> copy_output_until_stop [label=\"%d\"];\n", clava_dcg_global[52]);
   }
   if(clava_dcg_global[53] != 0) {
      fprintf(log_file_365, "	handle_compress -> isempty_RL [label=\"%d\"];\n", clava_dcg_global[53]);
   }
   if(clava_dcg_global[54] != 0) {
      fprintf(log_file_365, "	handle_compress -> prepare_new_block [label=\"%d\"];\n", clava_dcg_global[54]);
   }
   if(clava_dcg_global[55] != 0) {
      fprintf(log_file_365, "	handle_compress -> copy_input_until_stop [label=\"%d\"];\n", clava_dcg_global[55]);
   }
   if(clava_dcg_global[56] != 0) {
      fprintf(log_file_365, "	handle_compress -> flush_RL [label=\"%d\"];\n", clava_dcg_global[56]);
   }
   if(clava_dcg_global[57] != 0) {
      fprintf(log_file_365, "	handle_compress -> BZ2_compressBlock [label=\"%d\"];\n", clava_dcg_global[57]);
   }
   if(clava_dcg_global[58] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompress -> handle_compress [label=\"%d\"];\n", clava_dcg_global[58]);
   }
   if(clava_dcg_global[59] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompress -> isempty_RL [label=\"%d\"];\n", clava_dcg_global[59]);
   }
   if(clava_dcg_global[60] != 0) {
      fprintf(log_file_365, "	BZ2_bzCompressEnd -> strm [label=\"%d\"];\n", clava_dcg_global[60]);
   }
   if(clava_dcg_global[61] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompressInit -> bz_config_ok [label=\"%d\"];\n", clava_dcg_global[61]);
   }
   if(clava_dcg_global[62] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompressInit -> strm [label=\"%d\"];\n", clava_dcg_global[62]);
   }
   if(clava_dcg_global[63] != 0) {
      fprintf(log_file_365, "	unRLE_obuf_to_output_SMALL -> BZ2_indexIntoF [label=\"%d\"];\n", clava_dcg_global[63]);
   }
   if(clava_dcg_global[64] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompress -> unRLE_obuf_to_output_SMALL [label=\"%d\"];\n", clava_dcg_global[64]);
   }
   if(clava_dcg_global[65] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompress -> unRLE_obuf_to_output_FAST [label=\"%d\"];\n", clava_dcg_global[65]);
   }
   if(clava_dcg_global[66] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompress -> fprintf [label=\"%d\"];\n", clava_dcg_global[66]);
   }
   if(clava_dcg_global[67] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompress -> BZ2_decompress [label=\"%d\"];\n", clava_dcg_global[67]);
   }
   if(clava_dcg_global[68] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompress -> BZ2_bz__AssertH__fail [label=\"%d\"];\n", clava_dcg_global[68]);
   }
   if(clava_dcg_global[69] != 0) {
      fprintf(log_file_365, "	BZ2_bzDecompressEnd -> strm [label=\"%d\"];\n", clava_dcg_global[69]);
   }
   if(clava_dcg_global[70] != 0) {
      fprintf(log_file_365, "	myfeof -> fgetc [label=\"%d\"];\n", clava_dcg_global[70]);
   }
   if(clava_dcg_global[71] != 0) {
      fprintf(log_file_365, "	myfeof -> ungetc [label=\"%d\"];\n", clava_dcg_global[71]);
   }
   if(clava_dcg_global[72] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteOpen -> ferror [label=\"%d\"];\n", clava_dcg_global[72]);
   }
   if(clava_dcg_global[73] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteOpen -> malloc [label=\"%d\"];\n", clava_dcg_global[73]);
   }
   if(clava_dcg_global[74] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteOpen -> BZ2_bzCompressInit [label=\"%d\"];\n", clava_dcg_global[74]);
   }
   if(clava_dcg_global[75] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteOpen -> free [label=\"%d\"];\n", clava_dcg_global[75]);
   }
   if(clava_dcg_global[76] != 0) {
      fprintf(log_file_365, "	BZ2_bzWrite -> ferror [label=\"%d\"];\n", clava_dcg_global[76]);
   }
   if(clava_dcg_global[77] != 0) {
      fprintf(log_file_365, "	BZ2_bzWrite -> BZ2_bzCompress [label=\"%d\"];\n", clava_dcg_global[77]);
   }
   if(clava_dcg_global[78] != 0) {
      fprintf(log_file_365, "	BZ2_bzWrite -> fwrite [label=\"%d\"];\n", clava_dcg_global[78]);
   }
   if(clava_dcg_global[79] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose -> BZ2_bzWriteClose64 [label=\"%d\"];\n", clava_dcg_global[79]);
   }
   if(clava_dcg_global[80] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> ferror [label=\"%d\"];\n", clava_dcg_global[80]);
   }
   if(clava_dcg_global[81] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> BZ2_bzCompress [label=\"%d\"];\n", clava_dcg_global[81]);
   }
   if(clava_dcg_global[82] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> fwrite [label=\"%d\"];\n", clava_dcg_global[82]);
   }
   if(clava_dcg_global[83] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> fflush [label=\"%d\"];\n", clava_dcg_global[83]);
   }
   if(clava_dcg_global[84] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> BZ2_bzCompressEnd [label=\"%d\"];\n", clava_dcg_global[84]);
   }
   if(clava_dcg_global[85] != 0) {
      fprintf(log_file_365, "	BZ2_bzWriteClose64 -> free [label=\"%d\"];\n", clava_dcg_global[85]);
   }
   if(clava_dcg_global[86] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadOpen -> ferror [label=\"%d\"];\n", clava_dcg_global[86]);
   }
   if(clava_dcg_global[87] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadOpen -> malloc [label=\"%d\"];\n", clava_dcg_global[87]);
   }
   if(clava_dcg_global[88] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadOpen -> BZ2_bzDecompressInit [label=\"%d\"];\n", clava_dcg_global[88]);
   }
   if(clava_dcg_global[89] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadOpen -> free [label=\"%d\"];\n", clava_dcg_global[89]);
   }
   if(clava_dcg_global[90] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadClose -> BZ2_bzDecompressEnd [label=\"%d\"];\n", clava_dcg_global[90]);
   }
   if(clava_dcg_global[91] != 0) {
      fprintf(log_file_365, "	BZ2_bzReadClose -> free [label=\"%d\"];\n", clava_dcg_global[91]);
   }
   if(clava_dcg_global[92] != 0) {
      fprintf(log_file_365, "	BZ2_bzRead -> ferror [label=\"%d\"];\n", clava_dcg_global[92]);
   }
   if(clava_dcg_global[93] != 0) {
      fprintf(log_file_365, "	BZ2_bzRead -> myfeof [label=\"%d\"];\n", clava_dcg_global[93]);
   }
   if(clava_dcg_global[94] != 0) {
      fprintf(log_file_365, "	BZ2_bzRead -> fread [label=\"%d\"];\n", clava_dcg_global[94]);
   }
   if(clava_dcg_global[95] != 0) {
      fprintf(log_file_365, "	BZ2_bzRead -> BZ2_bzDecompress [label=\"%d\"];\n", clava_dcg_global[95]);
   }
   if(clava_dcg_global[96] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffCompress -> BZ2_bzCompressInit [label=\"%d\"];\n", clava_dcg_global[96]);
   }
   if(clava_dcg_global[97] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffCompress -> BZ2_bzCompress [label=\"%d\"];\n", clava_dcg_global[97]);
   }
   if(clava_dcg_global[98] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffCompress -> BZ2_bzCompressEnd [label=\"%d\"];\n", clava_dcg_global[98]);
   }
   if(clava_dcg_global[99] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffDecompress -> BZ2_bzDecompressInit [label=\"%d\"];\n", clava_dcg_global[99]);
   }
   if(clava_dcg_global[100] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffDecompress -> BZ2_bzDecompress [label=\"%d\"];\n", clava_dcg_global[100]);
   }
   if(clava_dcg_global[101] != 0) {
      fprintf(log_file_365, "	BZ2_bzBuffToBuffDecompress -> BZ2_bzDecompressEnd [label=\"%d\"];\n", clava_dcg_global[101]);
   }
   if(clava_dcg_global[102] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> __ctype_b_loc [label=\"%d\"];\n", clava_dcg_global[102]);
   }
   if(clava_dcg_global[103] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> strcat [label=\"%d\"];\n", clava_dcg_global[103]);
   }
   if(clava_dcg_global[104] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> strcmp [label=\"%d\"];\n", clava_dcg_global[104]);
   }
   if(clava_dcg_global[105] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> fopen [label=\"%d\"];\n", clava_dcg_global[105]);
   }
   if(clava_dcg_global[106] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> fdopen [label=\"%d\"];\n", clava_dcg_global[106]);
   }
   if(clava_dcg_global[107] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> BZ2_bzWriteOpen [label=\"%d\"];\n", clava_dcg_global[107]);
   }
   if(clava_dcg_global[108] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> BZ2_bzReadOpen [label=\"%d\"];\n", clava_dcg_global[108]);
   }
   if(clava_dcg_global[109] != 0) {
      fprintf(log_file_365, "	bzopen_or_bzdopen -> fclose [label=\"%d\"];\n", clava_dcg_global[109]);
   }
   if(clava_dcg_global[110] != 0) {
      fprintf(log_file_365, "	BZ2_bzopen -> bzopen_or_bzdopen [label=\"%d\"];\n", clava_dcg_global[110]);
   }
   if(clava_dcg_global[111] != 0) {
      fprintf(log_file_365, "	BZ2_bzdopen -> bzopen_or_bzdopen [label=\"%d\"];\n", clava_dcg_global[111]);
   }
   if(clava_dcg_global[112] != 0) {
      fprintf(log_file_365, "	BZ2_bzread -> BZ2_bzRead [label=\"%d\"];\n", clava_dcg_global[112]);
   }
   if(clava_dcg_global[113] != 0) {
      fprintf(log_file_365, "	BZ2_bzwrite -> BZ2_bzWrite [label=\"%d\"];\n", clava_dcg_global[113]);
   }
   if(clava_dcg_global[114] != 0) {
      fprintf(log_file_365, "	BZ2_bzclose -> BZ2_bzWriteClose [label=\"%d\"];\n", clava_dcg_global[114]);
   }
   if(clava_dcg_global[115] != 0) {
      fprintf(log_file_365, "	BZ2_bzclose -> BZ2_bzReadClose [label=\"%d\"];\n", clava_dcg_global[115]);
   }
   if(clava_dcg_global[116] != 0) {
      fprintf(log_file_365, "	BZ2_bzclose -> fclose [label=\"%d\"];\n", clava_dcg_global[116]);
   }
   if(clava_dcg_global[117] != 0) {
      fprintf(log_file_365, "	uInt64_toAscii -> uInt64_qrm10 [label=\"%d\"];\n", clava_dcg_global[117]);
   }
   if(clava_dcg_global[118] != 0) {
      fprintf(log_file_365, "	uInt64_toAscii -> uInt64_isZero [label=\"%d\"];\n", clava_dcg_global[118]);
   }
   if(clava_dcg_global[119] != 0) {
      fprintf(log_file_365, "	compressStream -> ferror [label=\"%d\"];\n", clava_dcg_global[119]);
   }
   if(clava_dcg_global[120] != 0) {
      fprintf(log_file_365, "	compressStream -> BZ2_bzWriteOpen [label=\"%d\"];\n", clava_dcg_global[120]);
   }
   if(clava_dcg_global[121] != 0) {
      fprintf(log_file_365, "	compressStream -> fprintf [label=\"%d\"];\n", clava_dcg_global[121]);
   }
   if(clava_dcg_global[122] != 0) {
      fprintf(log_file_365, "	compressStream -> myfeof [label=\"%d\"];\n", clava_dcg_global[122]);
   }
   if(clava_dcg_global[123] != 0) {
      fprintf(log_file_365, "	compressStream -> fread [label=\"%d\"];\n", clava_dcg_global[123]);
   }
   if(clava_dcg_global[124] != 0) {
      fprintf(log_file_365, "	compressStream -> BZ2_bzWrite [label=\"%d\"];\n", clava_dcg_global[124]);
   }
   if(clava_dcg_global[125] != 0) {
      fprintf(log_file_365, "	compressStream -> BZ2_bzWriteClose64 [label=\"%d\"];\n", clava_dcg_global[125]);
   }
   if(clava_dcg_global[126] != 0) {
      fprintf(log_file_365, "	compressStream -> fflush [label=\"%d\"];\n", clava_dcg_global[126]);
   }
   if(clava_dcg_global[127] != 0) {
      fprintf(log_file_365, "	compressStream -> fclose [label=\"%d\"];\n", clava_dcg_global[127]);
   }
   if(clava_dcg_global[128] != 0) {
      fprintf(log_file_365, "	compressStream -> uInt64_from_UInt32s [label=\"%d\"];\n", clava_dcg_global[128]);
   }
   if(clava_dcg_global[129] != 0) {
      fprintf(log_file_365, "	compressStream -> uInt64_to_double [label=\"%d\"];\n", clava_dcg_global[129]);
   }
   if(clava_dcg_global[130] != 0) {
      fprintf(log_file_365, "	compressStream -> uInt64_toAscii [label=\"%d\"];\n", clava_dcg_global[130]);
   }
   if(clava_dcg_global[131] != 0) {
      fprintf(log_file_365, "	compressStream -> configError [label=\"%d\"];\n", clava_dcg_global[131]);
   }
   if(clava_dcg_global[132] != 0) {
      fprintf(log_file_365, "	compressStream -> outOfMemory [label=\"%d\"];\n", clava_dcg_global[132]);
   }
   if(clava_dcg_global[133] != 0) {
      fprintf(log_file_365, "	compressStream -> ioError [label=\"%d\"];\n", clava_dcg_global[133]);
   }
   if(clava_dcg_global[134] != 0) {
      fprintf(log_file_365, "	compressStream -> panic [label=\"%d\"];\n", clava_dcg_global[134]);
   }
   if(clava_dcg_global[135] != 0) {
      fprintf(log_file_365, "	uncompressStream -> ferror [label=\"%d\"];\n", clava_dcg_global[135]);
   }
   if(clava_dcg_global[136] != 0) {
      fprintf(log_file_365, "	uncompressStream -> BZ2_bzReadOpen [label=\"%d\"];\n", clava_dcg_global[136]);
   }
   if(clava_dcg_global[137] != 0) {
      fprintf(log_file_365, "	uncompressStream -> BZ2_bzRead [label=\"%d\"];\n", clava_dcg_global[137]);
   }
   if(clava_dcg_global[138] != 0) {
      fprintf(log_file_365, "	uncompressStream -> fwrite [label=\"%d\"];\n", clava_dcg_global[138]);
   }
   if(clava_dcg_global[139] != 0) {
      fprintf(log_file_365, "	uncompressStream -> BZ2_bzReadGetUnused [label=\"%d\"];\n", clava_dcg_global[139]);
   }
   if(clava_dcg_global[140] != 0) {
      fprintf(log_file_365, "	uncompressStream -> panic [label=\"%d\"];\n", clava_dcg_global[140]);
   }
   if(clava_dcg_global[141] != 0) {
      fprintf(log_file_365, "	uncompressStream -> BZ2_bzReadClose [label=\"%d\"];\n", clava_dcg_global[141]);
   }
   if(clava_dcg_global[142] != 0) {
      fprintf(log_file_365, "	uncompressStream -> myfeof [label=\"%d\"];\n", clava_dcg_global[142]);
   }
   if(clava_dcg_global[143] != 0) {
      fprintf(log_file_365, "	uncompressStream -> fclose [label=\"%d\"];\n", clava_dcg_global[143]);
   }
   if(clava_dcg_global[144] != 0) {
      fprintf(log_file_365, "	uncompressStream -> fflush [label=\"%d\"];\n", clava_dcg_global[144]);
   }
   if(clava_dcg_global[145] != 0) {
      fprintf(log_file_365, "	uncompressStream -> fprintf [label=\"%d\"];\n", clava_dcg_global[145]);
   }
   if(clava_dcg_global[146] != 0) {
      fprintf(log_file_365, "	uncompressStream -> rewind [label=\"%d\"];\n", clava_dcg_global[146]);
   }
   if(clava_dcg_global[147] != 0) {
      fprintf(log_file_365, "	uncompressStream -> fread [label=\"%d\"];\n", clava_dcg_global[147]);
   }
   if(clava_dcg_global[148] != 0) {
      fprintf(log_file_365, "	uncompressStream -> configError [label=\"%d\"];\n", clava_dcg_global[148]);
   }
   if(clava_dcg_global[149] != 0) {
      fprintf(log_file_365, "	uncompressStream -> ioError [label=\"%d\"];\n", clava_dcg_global[149]);
   }
   if(clava_dcg_global[150] != 0) {
      fprintf(log_file_365, "	uncompressStream -> crcError [label=\"%d\"];\n", clava_dcg_global[150]);
   }
   if(clava_dcg_global[151] != 0) {
      fprintf(log_file_365, "	uncompressStream -> outOfMemory [label=\"%d\"];\n", clava_dcg_global[151]);
   }
   if(clava_dcg_global[152] != 0) {
      fprintf(log_file_365, "	uncompressStream -> compressedStreamEOF [label=\"%d\"];\n", clava_dcg_global[152]);
   }
   if(clava_dcg_global[153] != 0) {
      fprintf(log_file_365, "	testStream -> ferror [label=\"%d\"];\n", clava_dcg_global[153]);
   }
   if(clava_dcg_global[154] != 0) {
      fprintf(log_file_365, "	testStream -> BZ2_bzReadOpen [label=\"%d\"];\n", clava_dcg_global[154]);
   }
   if(clava_dcg_global[155] != 0) {
      fprintf(log_file_365, "	testStream -> BZ2_bzRead [label=\"%d\"];\n", clava_dcg_global[155]);
   }
   if(clava_dcg_global[156] != 0) {
      fprintf(log_file_365, "	testStream -> BZ2_bzReadGetUnused [label=\"%d\"];\n", clava_dcg_global[156]);
   }
   if(clava_dcg_global[157] != 0) {
      fprintf(log_file_365, "	testStream -> panic [label=\"%d\"];\n", clava_dcg_global[157]);
   }
   if(clava_dcg_global[158] != 0) {
      fprintf(log_file_365, "	testStream -> BZ2_bzReadClose [label=\"%d\"];\n", clava_dcg_global[158]);
   }
   if(clava_dcg_global[159] != 0) {
      fprintf(log_file_365, "	testStream -> myfeof [label=\"%d\"];\n", clava_dcg_global[159]);
   }
   if(clava_dcg_global[160] != 0) {
      fprintf(log_file_365, "	testStream -> fclose [label=\"%d\"];\n", clava_dcg_global[160]);
   }
   if(clava_dcg_global[161] != 0) {
      fprintf(log_file_365, "	testStream -> fprintf [label=\"%d\"];\n", clava_dcg_global[161]);
   }
   if(clava_dcg_global[162] != 0) {
      fprintf(log_file_365, "	testStream -> configError [label=\"%d\"];\n", clava_dcg_global[162]);
   }
   if(clava_dcg_global[163] != 0) {
      fprintf(log_file_365, "	testStream -> ioError [label=\"%d\"];\n", clava_dcg_global[163]);
   }
   if(clava_dcg_global[164] != 0) {
      fprintf(log_file_365, "	testStream -> outOfMemory [label=\"%d\"];\n", clava_dcg_global[164]);
   }
   if(clava_dcg_global[165] != 0) {
      fprintf(log_file_365, "	cadvise -> fprintf [label=\"%d\"];\n", clava_dcg_global[165]);
   }
   if(clava_dcg_global[166] != 0) {
      fprintf(log_file_365, "	showFileNames -> fprintf [label=\"%d\"];\n", clava_dcg_global[166]);
   }
   if(clava_dcg_global[167] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> stat [label=\"%d\"];\n", clava_dcg_global[167]);
   }
   if(clava_dcg_global[168] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> fprintf [label=\"%d\"];\n", clava_dcg_global[168]);
   }
   if(clava_dcg_global[169] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> fclose [label=\"%d\"];\n", clava_dcg_global[169]);
   }
   if(clava_dcg_global[170] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> remove [label=\"%d\"];\n", clava_dcg_global[170]);
   }
   if(clava_dcg_global[171] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> setExit [label=\"%d\"];\n", clava_dcg_global[171]);
   }
   if(clava_dcg_global[172] != 0) {
      fprintf(log_file_365, "	cleanUpAndFail -> exit [label=\"%d\"];\n", clava_dcg_global[172]);
   }
   if(clava_dcg_global[173] != 0) {
      fprintf(log_file_365, "	panic -> fprintf [label=\"%d\"];\n", clava_dcg_global[173]);
   }
   if(clava_dcg_global[174] != 0) {
      fprintf(log_file_365, "	panic -> showFileNames [label=\"%d\"];\n", clava_dcg_global[174]);
   }
   if(clava_dcg_global[175] != 0) {
      fprintf(log_file_365, "	panic -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[175]);
   }
   if(clava_dcg_global[176] != 0) {
      fprintf(log_file_365, "	crcError -> fprintf [label=\"%d\"];\n", clava_dcg_global[176]);
   }
   if(clava_dcg_global[177] != 0) {
      fprintf(log_file_365, "	crcError -> showFileNames [label=\"%d\"];\n", clava_dcg_global[177]);
   }
   if(clava_dcg_global[178] != 0) {
      fprintf(log_file_365, "	crcError -> cadvise [label=\"%d\"];\n", clava_dcg_global[178]);
   }
   if(clava_dcg_global[179] != 0) {
      fprintf(log_file_365, "	crcError -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[179]);
   }
   if(clava_dcg_global[180] != 0) {
      fprintf(log_file_365, "	compressedStreamEOF -> fprintf [label=\"%d\"];\n", clava_dcg_global[180]);
   }
   if(clava_dcg_global[181] != 0) {
      fprintf(log_file_365, "	compressedStreamEOF -> perror [label=\"%d\"];\n", clava_dcg_global[181]);
   }
   if(clava_dcg_global[182] != 0) {
      fprintf(log_file_365, "	compressedStreamEOF -> showFileNames [label=\"%d\"];\n", clava_dcg_global[182]);
   }
   if(clava_dcg_global[183] != 0) {
      fprintf(log_file_365, "	compressedStreamEOF -> cadvise [label=\"%d\"];\n", clava_dcg_global[183]);
   }
   if(clava_dcg_global[184] != 0) {
      fprintf(log_file_365, "	compressedStreamEOF -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[184]);
   }
   if(clava_dcg_global[185] != 0) {
      fprintf(log_file_365, "	ioError -> fprintf [label=\"%d\"];\n", clava_dcg_global[185]);
   }
   if(clava_dcg_global[186] != 0) {
      fprintf(log_file_365, "	ioError -> perror [label=\"%d\"];\n", clava_dcg_global[186]);
   }
   if(clava_dcg_global[187] != 0) {
      fprintf(log_file_365, "	ioError -> showFileNames [label=\"%d\"];\n", clava_dcg_global[187]);
   }
   if(clava_dcg_global[188] != 0) {
      fprintf(log_file_365, "	ioError -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[188]);
   }
   if(clava_dcg_global[189] != 0) {
      fprintf(log_file_365, "	mySignalCatcher -> fprintf [label=\"%d\"];\n", clava_dcg_global[189]);
   }
   if(clava_dcg_global[190] != 0) {
      fprintf(log_file_365, "	mySignalCatcher -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[190]);
   }
   if(clava_dcg_global[191] != 0) {
      fprintf(log_file_365, "	mySIGSEGVorSIGBUScatcher -> fprintf [label=\"%d\"];\n", clava_dcg_global[191]);
   }
   if(clava_dcg_global[192] != 0) {
      fprintf(log_file_365, "	mySIGSEGVorSIGBUScatcher -> showFileNames [label=\"%d\"];\n", clava_dcg_global[192]);
   }
   if(clava_dcg_global[193] != 0) {
      fprintf(log_file_365, "	mySIGSEGVorSIGBUScatcher -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[193]);
   }
   if(clava_dcg_global[194] != 0) {
      fprintf(log_file_365, "	mySIGSEGVorSIGBUScatcher -> cadvise [label=\"%d\"];\n", clava_dcg_global[194]);
   }
   if(clava_dcg_global[195] != 0) {
      fprintf(log_file_365, "	outOfMemory -> fprintf [label=\"%d\"];\n", clava_dcg_global[195]);
   }
   if(clava_dcg_global[196] != 0) {
      fprintf(log_file_365, "	outOfMemory -> showFileNames [label=\"%d\"];\n", clava_dcg_global[196]);
   }
   if(clava_dcg_global[197] != 0) {
      fprintf(log_file_365, "	outOfMemory -> cleanUpAndFail [label=\"%d\"];\n", clava_dcg_global[197]);
   }
   if(clava_dcg_global[198] != 0) {
      fprintf(log_file_365, "	configError -> fprintf [label=\"%d\"];\n", clava_dcg_global[198]);
   }
   if(clava_dcg_global[199] != 0) {
      fprintf(log_file_365, "	configError -> setExit [label=\"%d\"];\n", clava_dcg_global[199]);
   }
   if(clava_dcg_global[200] != 0) {
      fprintf(log_file_365, "	configError -> exit [label=\"%d\"];\n", clava_dcg_global[200]);
   }
   if(clava_dcg_global[201] != 0) {
      fprintf(log_file_365, "	pad -> strlen [label=\"%d\"];\n", clava_dcg_global[201]);
   }
   if(clava_dcg_global[202] != 0) {
      fprintf(log_file_365, "	pad -> fprintf [label=\"%d\"];\n", clava_dcg_global[202]);
   }
   if(clava_dcg_global[203] != 0) {
      fprintf(log_file_365, "	copyFileName -> strlen [label=\"%d\"];\n", clava_dcg_global[203]);
   }
   if(clava_dcg_global[204] != 0) {
      fprintf(log_file_365, "	copyFileName -> fprintf [label=\"%d\"];\n", clava_dcg_global[204]);
   }
   if(clava_dcg_global[205] != 0) {
      fprintf(log_file_365, "	copyFileName -> setExit [label=\"%d\"];\n", clava_dcg_global[205]);
   }
   if(clava_dcg_global[206] != 0) {
      fprintf(log_file_365, "	copyFileName -> exit [label=\"%d\"];\n", clava_dcg_global[206]);
   }
   if(clava_dcg_global[207] != 0) {
      fprintf(log_file_365, "	copyFileName -> strncpy [label=\"%d\"];\n", clava_dcg_global[207]);
   }
   if(clava_dcg_global[208] != 0) {
      fprintf(log_file_365, "	fileExists -> fopen [label=\"%d\"];\n", clava_dcg_global[208]);
   }
   if(clava_dcg_global[209] != 0) {
      fprintf(log_file_365, "	fileExists -> fclose [label=\"%d\"];\n", clava_dcg_global[209]);
   }
   if(clava_dcg_global[210] != 0) {
      fprintf(log_file_365, "	fopen_output_safely -> open [label=\"%d\"];\n", clava_dcg_global[210]);
   }
   if(clava_dcg_global[211] != 0) {
      fprintf(log_file_365, "	fopen_output_safely -> fdopen [label=\"%d\"];\n", clava_dcg_global[211]);
   }
   if(clava_dcg_global[212] != 0) {
      fprintf(log_file_365, "	fopen_output_safely -> close [label=\"%d\"];\n", clava_dcg_global[212]);
   }
   if(clava_dcg_global[213] != 0) {
      fprintf(log_file_365, "	notAStandardFile -> lstat [label=\"%d\"];\n", clava_dcg_global[213]);
   }
   if(clava_dcg_global[214] != 0) {
      fprintf(log_file_365, "	countHardLinks -> lstat [label=\"%d\"];\n", clava_dcg_global[214]);
   }
   if(clava_dcg_global[215] != 0) {
      fprintf(log_file_365, "	saveInputFileMetaInfo -> stat [label=\"%d\"];\n", clava_dcg_global[215]);
   }
   if(clava_dcg_global[216] != 0) {
      fprintf(log_file_365, "	saveInputFileMetaInfo -> ioError [label=\"%d\"];\n", clava_dcg_global[216]);
   }
   if(clava_dcg_global[217] != 0) {
      fprintf(log_file_365, "	applySavedMetaInfoToOutputFile -> chmod [label=\"%d\"];\n", clava_dcg_global[217]);
   }
   if(clava_dcg_global[218] != 0) {
      fprintf(log_file_365, "	applySavedMetaInfoToOutputFile -> ioError [label=\"%d\"];\n", clava_dcg_global[218]);
   }
   if(clava_dcg_global[219] != 0) {
      fprintf(log_file_365, "	applySavedMetaInfoToOutputFile -> utime [label=\"%d\"];\n", clava_dcg_global[219]);
   }
   if(clava_dcg_global[220] != 0) {
      fprintf(log_file_365, "	applySavedMetaInfoToOutputFile -> chown [label=\"%d\"];\n", clava_dcg_global[220]);
   }
   if(clava_dcg_global[221] != 0) {
      fprintf(log_file_365, "	hasSuffix -> strlen [label=\"%d\"];\n", clava_dcg_global[221]);
   }
   if(clava_dcg_global[222] != 0) {
      fprintf(log_file_365, "	hasSuffix -> strcmp [label=\"%d\"];\n", clava_dcg_global[222]);
   }
   if(clava_dcg_global[223] != 0) {
      fprintf(log_file_365, "	mapSuffix -> hasSuffix [label=\"%d\"];\n", clava_dcg_global[223]);
   }
   if(clava_dcg_global[224] != 0) {
      fprintf(log_file_365, "	mapSuffix -> strlen [label=\"%d\"];\n", clava_dcg_global[224]);
   }
   if(clava_dcg_global[225] != 0) {
      fprintf(log_file_365, "	mapSuffix -> strcat [label=\"%d\"];\n", clava_dcg_global[225]);
   }
   if(clava_dcg_global[226] != 0) {
      fprintf(log_file_365, "	compress -> panic [label=\"%d\"];\n", clava_dcg_global[226]);
   }
   if(clava_dcg_global[227] != 0) {
      fprintf(log_file_365, "	compress -> copyFileName [label=\"%d\"];\n", clava_dcg_global[227]);
   }
   if(clava_dcg_global[228] != 0) {
      fprintf(log_file_365, "	compress -> strcat [label=\"%d\"];\n", clava_dcg_global[228]);
   }
   if(clava_dcg_global[229] != 0) {
      fprintf(log_file_365, "	compress -> containsDubiousChars [label=\"%d\"];\n", clava_dcg_global[229]);
   }
   if(clava_dcg_global[230] != 0) {
      fprintf(log_file_365, "	compress -> fprintf [label=\"%d\"];\n", clava_dcg_global[230]);
   }
   if(clava_dcg_global[231] != 0) {
      fprintf(log_file_365, "	compress -> setExit [label=\"%d\"];\n", clava_dcg_global[231]);
   }
   if(clava_dcg_global[232] != 0) {
      fprintf(log_file_365, "	compress -> fileExists [label=\"%d\"];\n", clava_dcg_global[232]);
   }
   if(clava_dcg_global[233] != 0) {
      fprintf(log_file_365, "	compress -> strerror [label=\"%d\"];\n", clava_dcg_global[233]);
   }
   if(clava_dcg_global[234] != 0) {
      fprintf(log_file_365, "	compress -> __errno_location [label=\"%d\"];\n", clava_dcg_global[234]);
   }
   if(clava_dcg_global[235] != 0) {
      fprintf(log_file_365, "	compress -> hasSuffix [label=\"%d\"];\n", clava_dcg_global[235]);
   }
   if(clava_dcg_global[236] != 0) {
      fprintf(log_file_365, "	compress -> stat [label=\"%d\"];\n", clava_dcg_global[236]);
   }
   if(clava_dcg_global[237] != 0) {
      fprintf(log_file_365, "	compress -> notAStandardFile [label=\"%d\"];\n", clava_dcg_global[237]);
   }
   if(clava_dcg_global[238] != 0) {
      fprintf(log_file_365, "	compress -> remove [label=\"%d\"];\n", clava_dcg_global[238]);
   }
   if(clava_dcg_global[239] != 0) {
      fprintf(log_file_365, "	compress -> countHardLinks [label=\"%d\"];\n", clava_dcg_global[239]);
   }
   if(clava_dcg_global[240] != 0) {
      fprintf(log_file_365, "	compress -> saveInputFileMetaInfo [label=\"%d\"];\n", clava_dcg_global[240]);
   }
   if(clava_dcg_global[241] != 0) {
      fprintf(log_file_365, "	compress -> isatty [label=\"%d\"];\n", clava_dcg_global[241]);
   }
   if(clava_dcg_global[242] != 0) {
      fprintf(log_file_365, "	compress -> fileno [label=\"%d\"];\n", clava_dcg_global[242]);
   }
   if(clava_dcg_global[243] != 0) {
      fprintf(log_file_365, "	compress -> fopen [label=\"%d\"];\n", clava_dcg_global[243]);
   }
   if(clava_dcg_global[244] != 0) {
      fprintf(log_file_365, "	compress -> fclose [label=\"%d\"];\n", clava_dcg_global[244]);
   }
   if(clava_dcg_global[245] != 0) {
      fprintf(log_file_365, "	compress -> fopen_output_safely [label=\"%d\"];\n", clava_dcg_global[245]);
   }
   if(clava_dcg_global[246] != 0) {
      fprintf(log_file_365, "	compress -> pad [label=\"%d\"];\n", clava_dcg_global[246]);
   }
   if(clava_dcg_global[247] != 0) {
      fprintf(log_file_365, "	compress -> fflush [label=\"%d\"];\n", clava_dcg_global[247]);
   }
   if(clava_dcg_global[248] != 0) {
      fprintf(log_file_365, "	compress -> compressStream [label=\"%d\"];\n", clava_dcg_global[248]);
   }
   if(clava_dcg_global[249] != 0) {
      fprintf(log_file_365, "	compress -> applySavedMetaInfoToOutputFile [label=\"%d\"];\n", clava_dcg_global[249]);
   }
   if(clava_dcg_global[250] != 0) {
      fprintf(log_file_365, "	compress -> ioError [label=\"%d\"];\n", clava_dcg_global[250]);
   }
   if(clava_dcg_global[251] != 0) {
      fprintf(log_file_365, "	uncompress -> panic [label=\"%d\"];\n", clava_dcg_global[251]);
   }
   if(clava_dcg_global[252] != 0) {
      fprintf(log_file_365, "	uncompress -> copyFileName [label=\"%d\"];\n", clava_dcg_global[252]);
   }
   if(clava_dcg_global[253] != 0) {
      fprintf(log_file_365, "	uncompress -> mapSuffix [label=\"%d\"];\n", clava_dcg_global[253]);
   }
   if(clava_dcg_global[254] != 0) {
      fprintf(log_file_365, "	uncompress -> strcat [label=\"%d\"];\n", clava_dcg_global[254]);
   }
   if(clava_dcg_global[255] != 0) {
      fprintf(log_file_365, "	uncompress -> containsDubiousChars [label=\"%d\"];\n", clava_dcg_global[255]);
   }
   if(clava_dcg_global[256] != 0) {
      fprintf(log_file_365, "	uncompress -> fprintf [label=\"%d\"];\n", clava_dcg_global[256]);
   }
   if(clava_dcg_global[257] != 0) {
      fprintf(log_file_365, "	uncompress -> setExit [label=\"%d\"];\n", clava_dcg_global[257]);
   }
   if(clava_dcg_global[258] != 0) {
      fprintf(log_file_365, "	uncompress -> fileExists [label=\"%d\"];\n", clava_dcg_global[258]);
   }
   if(clava_dcg_global[259] != 0) {
      fprintf(log_file_365, "	uncompress -> strerror [label=\"%d\"];\n", clava_dcg_global[259]);
   }
   if(clava_dcg_global[260] != 0) {
      fprintf(log_file_365, "	uncompress -> __errno_location [label=\"%d\"];\n", clava_dcg_global[260]);
   }
   if(clava_dcg_global[261] != 0) {
      fprintf(log_file_365, "	uncompress -> stat [label=\"%d\"];\n", clava_dcg_global[261]);
   }
   if(clava_dcg_global[262] != 0) {
      fprintf(log_file_365, "	uncompress -> notAStandardFile [label=\"%d\"];\n", clava_dcg_global[262]);
   }
   if(clava_dcg_global[263] != 0) {
      fprintf(log_file_365, "	uncompress -> remove [label=\"%d\"];\n", clava_dcg_global[263]);
   }
   if(clava_dcg_global[264] != 0) {
      fprintf(log_file_365, "	uncompress -> countHardLinks [label=\"%d\"];\n", clava_dcg_global[264]);
   }
   if(clava_dcg_global[265] != 0) {
      fprintf(log_file_365, "	uncompress -> saveInputFileMetaInfo [label=\"%d\"];\n", clava_dcg_global[265]);
   }
   if(clava_dcg_global[266] != 0) {
      fprintf(log_file_365, "	uncompress -> isatty [label=\"%d\"];\n", clava_dcg_global[266]);
   }
   if(clava_dcg_global[267] != 0) {
      fprintf(log_file_365, "	uncompress -> fileno [label=\"%d\"];\n", clava_dcg_global[267]);
   }
   if(clava_dcg_global[268] != 0) {
      fprintf(log_file_365, "	uncompress -> fopen [label=\"%d\"];\n", clava_dcg_global[268]);
   }
   if(clava_dcg_global[269] != 0) {
      fprintf(log_file_365, "	uncompress -> fclose [label=\"%d\"];\n", clava_dcg_global[269]);
   }
   if(clava_dcg_global[270] != 0) {
      fprintf(log_file_365, "	uncompress -> fopen_output_safely [label=\"%d\"];\n", clava_dcg_global[270]);
   }
   if(clava_dcg_global[271] != 0) {
      fprintf(log_file_365, "	uncompress -> pad [label=\"%d\"];\n", clava_dcg_global[271]);
   }
   if(clava_dcg_global[272] != 0) {
      fprintf(log_file_365, "	uncompress -> fflush [label=\"%d\"];\n", clava_dcg_global[272]);
   }
   if(clava_dcg_global[273] != 0) {
      fprintf(log_file_365, "	uncompress -> uncompressStream [label=\"%d\"];\n", clava_dcg_global[273]);
   }
   if(clava_dcg_global[274] != 0) {
      fprintf(log_file_365, "	uncompress -> applySavedMetaInfoToOutputFile [label=\"%d\"];\n", clava_dcg_global[274]);
   }
   if(clava_dcg_global[275] != 0) {
      fprintf(log_file_365, "	uncompress -> ioError [label=\"%d\"];\n", clava_dcg_global[275]);
   }
   if(clava_dcg_global[276] != 0) {
      fprintf(log_file_365, "	testf -> panic [label=\"%d\"];\n", clava_dcg_global[276]);
   }
   if(clava_dcg_global[277] != 0) {
      fprintf(log_file_365, "	testf -> copyFileName [label=\"%d\"];\n", clava_dcg_global[277]);
   }
   if(clava_dcg_global[278] != 0) {
      fprintf(log_file_365, "	testf -> containsDubiousChars [label=\"%d\"];\n", clava_dcg_global[278]);
   }
   if(clava_dcg_global[279] != 0) {
      fprintf(log_file_365, "	testf -> fprintf [label=\"%d\"];\n", clava_dcg_global[279]);
   }
   if(clava_dcg_global[280] != 0) {
      fprintf(log_file_365, "	testf -> setExit [label=\"%d\"];\n", clava_dcg_global[280]);
   }
   if(clava_dcg_global[281] != 0) {
      fprintf(log_file_365, "	testf -> fileExists [label=\"%d\"];\n", clava_dcg_global[281]);
   }
   if(clava_dcg_global[282] != 0) {
      fprintf(log_file_365, "	testf -> strerror [label=\"%d\"];\n", clava_dcg_global[282]);
   }
   if(clava_dcg_global[283] != 0) {
      fprintf(log_file_365, "	testf -> __errno_location [label=\"%d\"];\n", clava_dcg_global[283]);
   }
   if(clava_dcg_global[284] != 0) {
      fprintf(log_file_365, "	testf -> stat [label=\"%d\"];\n", clava_dcg_global[284]);
   }
   if(clava_dcg_global[285] != 0) {
      fprintf(log_file_365, "	testf -> isatty [label=\"%d\"];\n", clava_dcg_global[285]);
   }
   if(clava_dcg_global[286] != 0) {
      fprintf(log_file_365, "	testf -> fileno [label=\"%d\"];\n", clava_dcg_global[286]);
   }
   if(clava_dcg_global[287] != 0) {
      fprintf(log_file_365, "	testf -> fopen [label=\"%d\"];\n", clava_dcg_global[287]);
   }
   if(clava_dcg_global[288] != 0) {
      fprintf(log_file_365, "	testf -> pad [label=\"%d\"];\n", clava_dcg_global[288]);
   }
   if(clava_dcg_global[289] != 0) {
      fprintf(log_file_365, "	testf -> fflush [label=\"%d\"];\n", clava_dcg_global[289]);
   }
   if(clava_dcg_global[290] != 0) {
      fprintf(log_file_365, "	testf -> testStream [label=\"%d\"];\n", clava_dcg_global[290]);
   }
   if(clava_dcg_global[291] != 0) {
      fprintf(log_file_365, "	license -> fprintf [label=\"%d\"];\n", clava_dcg_global[291]);
   }
   if(clava_dcg_global[292] != 0) {
      fprintf(log_file_365, "	license -> BZ2_bzlibVersion [label=\"%d\"];\n", clava_dcg_global[292]);
   }
   if(clava_dcg_global[293] != 0) {
      fprintf(log_file_365, "	usage -> fprintf [label=\"%d\"];\n", clava_dcg_global[293]);
   }
   if(clava_dcg_global[294] != 0) {
      fprintf(log_file_365, "	usage -> BZ2_bzlibVersion [label=\"%d\"];\n", clava_dcg_global[294]);
   }
   if(clava_dcg_global[295] != 0) {
      fprintf(log_file_365, "	redundant -> fprintf [label=\"%d\"];\n", clava_dcg_global[295]);
   }
   if(clava_dcg_global[296] != 0) {
      fprintf(log_file_365, "	myMalloc -> malloc [label=\"%d\"];\n", clava_dcg_global[296]);
   }
   if(clava_dcg_global[297] != 0) {
      fprintf(log_file_365, "	myMalloc -> outOfMemory [label=\"%d\"];\n", clava_dcg_global[297]);
   }
   if(clava_dcg_global[298] != 0) {
      fprintf(log_file_365, "	mkCell -> myMalloc [label=\"%d\"];\n", clava_dcg_global[298]);
   }
   if(clava_dcg_global[299] != 0) {
      fprintf(log_file_365, "	snocString -> mkCell [label=\"%d\"];\n", clava_dcg_global[299]);
   }
   if(clava_dcg_global[300] != 0) {
      fprintf(log_file_365, "	snocString -> myMalloc [label=\"%d\"];\n", clava_dcg_global[300]);
   }
   if(clava_dcg_global[301] != 0) {
      fprintf(log_file_365, "	snocString -> strlen [label=\"%d\"];\n", clava_dcg_global[301]);
   }
   if(clava_dcg_global[302] != 0) {
      fprintf(log_file_365, "	snocString -> strcpy [label=\"%d\"];\n", clava_dcg_global[302]);
   }
   if(clava_dcg_global[303] != 0) {
      fprintf(log_file_365, "	snocString -> snocString [label=\"%d\"];\n", clava_dcg_global[303]);
   }
   if(clava_dcg_global[304] != 0) {
      fprintf(log_file_365, "	addFlagsFromEnvVar -> getenv [label=\"%d\"];\n", clava_dcg_global[304]);
   }
   if(clava_dcg_global[305] != 0) {
      fprintf(log_file_365, "	addFlagsFromEnvVar -> __ctype_b_loc [label=\"%d\"];\n", clava_dcg_global[305]);
   }
   if(clava_dcg_global[306] != 0) {
      fprintf(log_file_365, "	addFlagsFromEnvVar -> snocString [label=\"%d\"];\n", clava_dcg_global[306]);
   }
   if(clava_dcg_global[307] != 0) {
      fprintf(log_file_365, "	main -> configError [label=\"%d\"];\n", clava_dcg_global[307]);
   }
   if(clava_dcg_global[308] != 0) {
      fprintf(log_file_365, "	main -> signal [label=\"%d\"];\n", clava_dcg_global[308]);
   }
   if(clava_dcg_global[309] != 0) {
      fprintf(log_file_365, "	main -> copyFileName [label=\"%d\"];\n", clava_dcg_global[309]);
   }
   if(clava_dcg_global[310] != 0) {
      fprintf(log_file_365, "	main -> addFlagsFromEnvVar [label=\"%d\"];\n", clava_dcg_global[310]);
   }
   if(clava_dcg_global[311] != 0) {
      fprintf(log_file_365, "	main -> snocString [label=\"%d\"];\n", clava_dcg_global[311]);
   }
   if(clava_dcg_global[312] != 0) {
      fprintf(log_file_365, "	main -> strcmp [label=\"%d\"];\n", clava_dcg_global[312]);
   }
   if(clava_dcg_global[313] != 0) {
      fprintf(log_file_365, "	main -> strlen [label=\"%d\"];\n", clava_dcg_global[313]);
   }
   if(clava_dcg_global[314] != 0) {
      fprintf(log_file_365, "	main -> strstr [label=\"%d\"];\n", clava_dcg_global[314]);
   }
   if(clava_dcg_global[315] != 0) {
      fprintf(log_file_365, "	main -> license [label=\"%d\"];\n", clava_dcg_global[315]);
   }
   if(clava_dcg_global[316] != 0) {
      fprintf(log_file_365, "	main -> usage [label=\"%d\"];\n", clava_dcg_global[316]);
   }
   if(clava_dcg_global[317] != 0) {
      fprintf(log_file_365, "	main -> exit [label=\"%d\"];\n", clava_dcg_global[317]);
   }
   if(clava_dcg_global[318] != 0) {
      fprintf(log_file_365, "	main -> fprintf [label=\"%d\"];\n", clava_dcg_global[318]);
   }
   if(clava_dcg_global[319] != 0) {
      fprintf(log_file_365, "	main -> redundant [label=\"%d\"];\n", clava_dcg_global[319]);
   }
   if(clava_dcg_global[320] != 0) {
      fprintf(log_file_365, "	main -> strncmp [label=\"%d\"];\n", clava_dcg_global[320]);
   }
   if(clava_dcg_global[321] != 0) {
      fprintf(log_file_365, "	main -> compress [label=\"%d\"];\n", clava_dcg_global[321]);
   }
   if(clava_dcg_global[322] != 0) {
      fprintf(log_file_365, "	main -> uncompress [label=\"%d\"];\n", clava_dcg_global[322]);
   }
   if(clava_dcg_global[323] != 0) {
      fprintf(log_file_365, "	main -> setExit [label=\"%d\"];\n", clava_dcg_global[323]);
   }
   if(clava_dcg_global[324] != 0) {
      fprintf(log_file_365, "	main -> testf [label=\"%d\"];\n", clava_dcg_global[324]);
   }
   if(clava_dcg_global[325] != 0) {
      fprintf(log_file_365, "	main -> free [label=\"%d\"];\n", clava_dcg_global[325]);
   }
   fprintf(log_file_365, "}\n");
   fclose(log_file_365);
}

IntNative main(IntNative argc, Char *argv[]) {
   atexit(clava_call_graph);
   Int32 i, j;
   Char *tmp;
   Cell *argList;
   Cell *aa;
   Bool decode;
   /*-- Be really really really paranoid :-) --*/
   if(sizeof(Int32) != 4 || sizeof(UInt32) != 4 || sizeof(Int16) != 2 || sizeof(UInt16) != 2 || sizeof(Char) != 1 || sizeof(UChar) != 1) {
      clava_dcg_global[ 307 ]++;
      configError();
   }
   /*-- Initialise --*/
   outputHandleJustInCase = ((void *) 0);
   smallMode = ((Bool) 0);
   keepInputFiles = ((Bool) 0);
   forceOverwrite = ((Bool) 0);
   noisy = ((Bool) 1);
   verbosity = 0;
   blockSize100k = 9;
   testFailsExist = ((Bool) 0);
   unzFailsExist = ((Bool) 0);
   numFileNames = 0;
   numFilesProcessed = 0;
   workFactor = 30;
   deleteOutputOnInterrupt = ((Bool) 0);
   exitValue = 0;
   i = j = 0;
   /*avoid bogus warning from egcs-1.1.X*/
   /*-- Set up signal handlers for mem access errors --*/
   clava_dcg_global[ 308 ]++;
   signal(11, mySIGSEGVorSIGBUScatcher);
   clava_dcg_global[ 308 ]++;
   signal(7, mySIGSEGVorSIGBUScatcher);
   clava_dcg_global[ 309 ]++;
   copyFileName(inName, "(none)");
   clava_dcg_global[ 309 ]++;
   copyFileName(outName, "(none)");
   clava_dcg_global[ 309 ]++;
   copyFileName(progNameReally, argv[0]);
   progName = &progNameReally[0];
   for(tmp = &progNameReally[0]; *tmp != '\0'; tmp++) if(*tmp == '/') progName = tmp + 1;
   /*-- Copy flags from env var BZIP2, and
   expand filename wildcards in arg list.
   --*/
   argList = ((void *) 0);
   clava_dcg_global[ 310 ]++;
   addFlagsFromEnvVar(&argList, "BZIP2");
   clava_dcg_global[ 310 ]++;
   addFlagsFromEnvVar(&argList, "BZIP");
   for(i = 1; i <= argc - 1; i++) {
      clava_dcg_global[ 311 ]++;
      argList = snocString((argList), (argv[i]));
   }
   /*-- Find the length of the longest filename --*/
   longestFileName = 7;
   numFileNames = 0;
   decode = ((Bool) 1);
   for(aa = argList; aa != ((void *) 0); aa = aa->link) {
      clava_dcg_global[ 312 ]++;
      if((strcmp(aa->name, ("--")) == 0)) {
         decode = ((Bool) 0);
         continue;
      }
      if(aa->name[0] == '-' && decode) continue;
      numFileNames++;
      clava_dcg_global[ 313 ]++;
      if(longestFileName < (Int32) strlen(aa->name)) {
         clava_dcg_global[ 313 ]++;
         longestFileName = (Int32) strlen(aa->name);
      }
   }
   /*-- Determine source modes; flag handling may change this too. --*/
   if(numFileNames == 0) srcMode = 1;
   else srcMode = 3;
   /*-- Determine what to do (compress/uncompress/test/cat). --*/
   /*-- Note that subsequent flag handling may change this. --*/
   opMode = 1;
   clava_dcg_global[ 314 ]++;
   clava_dcg_global[ 314 ]++;
   if((strstr(progName, "unzip") != 0) || (strstr(progName, "UNZIP") != 0)) opMode = 2;
   clava_dcg_global[ 314 ]++;
   clava_dcg_global[ 314 ]++;
   clava_dcg_global[ 314 ]++;
   clava_dcg_global[ 314 ]++;
   if((strstr(progName, "z2cat") != 0) || (strstr(progName, "Z2CAT") != 0) || (strstr(progName, "zcat") != 0) || (strstr(progName, "ZCAT") != 0)) {
      opMode = 2;
      srcMode = (numFileNames == 0) ? 1 : 2;
   }
   /*-- Look at the flags. --*/
   for(aa = argList; aa != ((void *) 0); aa = aa->link) {
      clava_dcg_global[ 312 ]++;
      if((strcmp(aa->name, ("--")) == 0)) break;
      if(aa->name[0] == '-' && aa->name[1] != '-') {
         for(j = 1; aa->name[j] != '\0'; j++) {
            switch (aa->name[j]) {
               case 'c':
               srcMode = 2;
               break;
               case 'd':
               opMode = 2;
               break;
               case 'z':
               opMode = 1;
               break;
               case 'f':
               forceOverwrite = ((Bool) 1);
               break;
               case 't':
               opMode = 3;
               break;
               case 'k':
               keepInputFiles = ((Bool) 1);
               break;
               case 's':
               smallMode = ((Bool) 1);
               break;
               case 'q':
               noisy = ((Bool) 0);
               break;
               case '1':
               blockSize100k = 1;
               break;
               case '2':
               blockSize100k = 2;
               break;
               case '3':
               blockSize100k = 3;
               break;
               case '4':
               blockSize100k = 4;
               break;
               case '5':
               blockSize100k = 5;
               break;
               case '6':
               blockSize100k = 6;
               break;
               case '7':
               blockSize100k = 7;
               break;
               case '8':
               blockSize100k = 8;
               break;
               case '9':
               blockSize100k = 9;
               break;
               case 'V':
               case 'L':
               clava_dcg_global[ 315 ]++;
               license();
               break;
               case 'v':
               verbosity++;
               break;
               case 'h':
               clava_dcg_global[ 316 ]++;
               usage(progName);
               clava_dcg_global[ 317 ]++;
               exit(0);
               break;
               default:
               clava_dcg_global[ 318 ]++;
               fprintf(stderr, "%s: Bad flag `%s'\n", progName, aa->name);
               clava_dcg_global[ 316 ]++;
               usage(progName);
               clava_dcg_global[ 317 ]++;
               exit(1);
               break;
            }
         }
      }
   }
   /*-- And again ... --*/
   for(aa = argList; aa != ((void *) 0); aa = aa->link) {
      clava_dcg_global[ 312 ]++;
      if((strcmp(aa->name, ("--")) == 0)) break;
      clava_dcg_global[ 312 ]++;
      if((strcmp(aa->name, ("--stdout")) == 0)) srcMode = 2;
      else {
         clava_dcg_global[ 312 ]++;
         if((strcmp(aa->name, ("--decompress")) == 0)) opMode = 2;
         else {
            clava_dcg_global[ 312 ]++;
            if((strcmp(aa->name, ("--compress")) == 0)) opMode = 1;
            else {
               clava_dcg_global[ 312 ]++;
               if((strcmp(aa->name, ("--force")) == 0)) forceOverwrite = ((Bool) 1);
               else {
                  clava_dcg_global[ 312 ]++;
                  if((strcmp(aa->name, ("--test")) == 0)) opMode = 3;
                  else {
                     clava_dcg_global[ 312 ]++;
                     if((strcmp(aa->name, ("--keep")) == 0)) keepInputFiles = ((Bool) 1);
                     else {
                        clava_dcg_global[ 312 ]++;
                        if((strcmp(aa->name, ("--small")) == 0)) smallMode = ((Bool) 1);
                        else {
                           clava_dcg_global[ 312 ]++;
                           if((strcmp(aa->name, ("--quiet")) == 0)) noisy = ((Bool) 0);
                           else {
                              clava_dcg_global[ 312 ]++;
                              if((strcmp(aa->name, ("--version")) == 0)) {
                                 clava_dcg_global[ 315 ]++;
                                 license();
                              }
                              else {
                                 clava_dcg_global[ 312 ]++;
                                 if((strcmp(aa->name, ("--license")) == 0)) {
                                    clava_dcg_global[ 315 ]++;
                                    license();
                                 }
                                 else {
                                    clava_dcg_global[ 312 ]++;
                                    if((strcmp(aa->name, ("--exponential")) == 0)) workFactor = 1;
                                    else {
                                       clava_dcg_global[ 312 ]++;
                                       if((strcmp(aa->name, ("--repetitive-best")) == 0)) {
                                          clava_dcg_global[ 319 ]++;
                                          redundant(aa->name);
                                       }
                                       else {
                                          clava_dcg_global[ 312 ]++;
                                          if((strcmp(aa->name, ("--repetitive-fast")) == 0)) {
                                             clava_dcg_global[ 319 ]++;
                                             redundant(aa->name);
                                          }
                                          else {
                                             clava_dcg_global[ 312 ]++;
                                             if((strcmp(aa->name, ("--fast")) == 0)) blockSize100k = 1;
                                             else {
                                                clava_dcg_global[ 312 ]++;
                                                if((strcmp(aa->name, ("--best")) == 0)) blockSize100k = 9;
                                                else {
                                                   clava_dcg_global[ 312 ]++;
                                                   if((strcmp(aa->name, ("--verbose")) == 0)) verbosity++;
                                                   else {
                                                      clava_dcg_global[ 312 ]++;
                                                      if((strcmp(aa->name, ("--help")) == 0)) {
                                                         clava_dcg_global[ 316 ]++;
                                                         usage(progName);
                                                         clava_dcg_global[ 317 ]++;
                                                         exit(0);
                                                      }
                                                      else {
                                                         clava_dcg_global[ 320 ]++;
                                                         if(strncmp(aa->name, "--", 2) == 0) {
                                                            clava_dcg_global[ 318 ]++;
                                                            fprintf(stderr, "%s: Bad flag `%s'\n", progName, aa->name);
                                                            clava_dcg_global[ 316 ]++;
                                                            usage(progName);
                                                            clava_dcg_global[ 317 ]++;
                                                            exit(1);
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if(verbosity > 4) verbosity = 4;
   if(opMode == 1 && smallMode && blockSize100k > 2) blockSize100k = 2;
   if(opMode == 3 && srcMode == 2) {
      clava_dcg_global[ 318 ]++;
      fprintf(stderr, "%s: -c and -t cannot be used together.\n", progName);
      clava_dcg_global[ 317 ]++;
      exit(1);
   }
   if(srcMode == 2 && numFileNames == 0) srcMode = 1;
   if(opMode != 1) blockSize100k = 0;
   if(srcMode == 3) {
      clava_dcg_global[ 308 ]++;
      signal(2, mySignalCatcher);
      clava_dcg_global[ 308 ]++;
      signal(15, mySignalCatcher);
      clava_dcg_global[ 308 ]++;
      signal(1, mySignalCatcher);
   }
   if(opMode == 1) {
      if(srcMode == 1) {
         clava_dcg_global[ 321 ]++;
         compress(((void *) 0));
      }
      else {
         decode = ((Bool) 1);
         for(aa = argList; aa != ((void *) 0); aa = aa->link) {
            clava_dcg_global[ 312 ]++;
            if((strcmp(aa->name, ("--")) == 0)) {
               decode = ((Bool) 0);
               continue;
            }
            if(aa->name[0] == '-' && decode) continue;
            numFilesProcessed++;
            clava_dcg_global[ 321 ]++;
            compress(aa->name);
         }
      }
   }
   else if(opMode == 2) {
      unzFailsExist = ((Bool) 0);
      if(srcMode == 1) {
         clava_dcg_global[ 322 ]++;
         uncompress(((void *) 0));
      }
      else {
         decode = ((Bool) 1);
         for(aa = argList; aa != ((void *) 0); aa = aa->link) {
            clava_dcg_global[ 312 ]++;
            if((strcmp(aa->name, ("--")) == 0)) {
               decode = ((Bool) 0);
               continue;
            }
            if(aa->name[0] == '-' && decode) continue;
            numFilesProcessed++;
            clava_dcg_global[ 322 ]++;
            uncompress(aa->name);
         }
      }
      if(unzFailsExist) {
         clava_dcg_global[ 323 ]++;
         setExit(2);
         clava_dcg_global[ 317 ]++;
         exit(exitValue);
      }
   }
   else {
      testFailsExist = ((Bool) 0);
      if(srcMode == 1) {
         clava_dcg_global[ 324 ]++;
         testf(((void *) 0));
      }
      else {
         decode = ((Bool) 1);
         for(aa = argList; aa != ((void *) 0); aa = aa->link) {
            clava_dcg_global[ 312 ]++;
            if((strcmp(aa->name, ("--")) == 0)) {
               decode = ((Bool) 0);
               continue;
            }
            if(aa->name[0] == '-' && decode) continue;
            numFilesProcessed++;
            clava_dcg_global[ 324 ]++;
            testf(aa->name);
         }
      }
      if(testFailsExist && noisy) {
         clava_dcg_global[ 318 ]++;
         fprintf(stderr, "\nYou can use the `bzip2recover' program to attempt to recover\ndata from undamaged sections of corrupted files.\n\n");
         clava_dcg_global[ 323 ]++;
         setExit(2);
         clava_dcg_global[ 317 ]++;
         exit(exitValue);
      }
   }
   /*Free the argument list memory to mollify leak detectors
   (eg) Purify, Checker.  Serves no other useful purpose.
   */
   aa = argList;
   while(aa != ((void *) 0)) {
      Cell *aa2 = aa->link;
      if(aa->name != ((void *) 0)) {
         clava_dcg_global[ 325 ]++;
         free(aa->name);
      }
      clava_dcg_global[ 325 ]++;
      free(aa);
      aa = aa2;
   }
   
   return exitValue;
}

/*-----------------------------------------------------------*/
/*--- end                                         bzip2.c ---*/
/*-----------------------------------------------------------*/