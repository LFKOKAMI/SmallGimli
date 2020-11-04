#pragma once
#ifndef GIMLI_H_
#define GIMLI_H_

#include <vector>
using namespace std;

typedef unsigned int uint32;

//used for finding valid capacity part
struct CapHalf {
	uint32 in[4];
	uint32 row[3];
};

const uint32 E[32] = {
	0x1,0x2,0x4,0x8,
	0x10,0x20,0x40,0x80,
	0x100,0x200,0x400,0x800,
	0x1000,0x2000,0x4000,0x8000,
	0x10000,0x20000,0x40000,0x80000,
	0x100000,0x200000,0x400000,0x800000,
	0x1000000,0x2000000,0x4000000,0x8000000,
	0x10000000,0x20000000,0x40000000,0x80000000
};

class Gimli {
private:
	int previousINT32;
	uint32 constant[6];
	vector<uint32> p3;
public:
	Gimli();
	~Gimli();
	uint32 getRand32();
	uint32 getRand16();
	uint32 getRand8();
	bool BIT(uint32 x, int n);
	void toBit(uint32 s, bool b[], int size);
	uint32 toUint32(bool b[], int size);
	uint32 RL(uint32 word, int n,int wordSize);
	uint32 L(uint32 word, int n,int wordSize);
	void SP(uint32 ix, uint32 iy, uint32 iz, uint32 &ox, uint32 &oy, uint32 &oz, int wordSize);
	void SPInver(uint32 input[], int wordSize);

	void bigSwap(uint32 state[]);
	void smallSwap(uint32 state[]);
	void SPColumn(uint32 state[], int size, int wordSize);
	void inverseSPColumn(uint32 state[], int size, int wordSize);

	void check5RPreimageAttack_FindCapacityPart();
	bool fillLeftWords_FindCap(CapHalf* list, int index, uint32 state[],uint32 hash[][4],uint32 s5filled[]);
	bool property3(uint32 ix, uint32 &iy, uint32 &iz,uint32 &ox,uint32 oy,uint32 oz);

	bool check5RPreimageAttack_MatchCapacityPart();
	void constructMappingTable(vector<vector<CapHalf> >& table, uint32 constant[], int num);
	bool property4(uint32 iy, uint32 iz, uint32 oy, uint32 oz, uint32& ox);
	void permutation5R(uint32 s[], uint32 constant[], int wordSize);
	void inversePermutation5R(uint32 s[], uint32 constant[],int wordSize);

	bool check9RPreimageAttackXOF_FindRate();
	bool fourRoundsOneColumn(uint32 y, uint32 z, uint32 outputX, CapHalf* table, int tableSize, CapHalf* capHalf, int size,int &finalSize, bool isSearch,uint32 &r0,uint32 &r1);
	void twoRoundsSP(uint32 outputX, CapHalf* table, int &size);

	uint32 checkFullRoundDis();
	bool zeroInternalDiffAttack18R();
	void permutation8HalfR(uint32 s[], uint32 constant[], int wordSize);

	void outputState(uint32 state[]);
	void permutation(uint32 state[], uint32 constant[],int wordSize);
};

#endif // !GIMLI_H_
