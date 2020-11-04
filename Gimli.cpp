#include "Gimli.h"
#include <ctime>
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace std;

void copyState(uint32 src[], uint32 des[]) {
	for (int i = 0; i < 12; i++) {
		des[i] = src[i];
	}
}

bool compare(CapHalf c0, CapHalf c1) {
	if (c0.row[2] > c1.row[2]) {
		return true;
	}
	return false;
}

int binary_search_find_index(CapHalf* v, uint32 size,CapHalf data) {
	auto it = std::lower_bound(v, v+size, data, compare);
	std::size_t index = std::distance(v, it);
	return index;
}

Gimli::Gimli() {
	//wordSize = 32;
	previousINT32 = 0;

	for (int i = 24; i > 0; i--) {
		if (i % 4 == 0) {
			constant[6 - i / 4] = 0x9e377900 ^ i;
		}
	}

	srand(time(NULL));
}

Gimli::~Gimli() {

}

void Gimli::permutation(uint32 state[], uint32 constant[],int wordSize) {
	for (int i = 24; i > 0; i--) {
		SPColumn(state, 12,wordSize);
		if (i % 4 == 0) {//small swap
			swap(state[0], state[3]);
			swap(state[6], state[9]);
		}
		else if (i % 2 == 0) {//Big swap
			swap(state[0], state[6]);
			swap(state[3], state[9]);
		}

		if (i % 4 == 0) {
			state[0] = state[0] ^ constant[6-i/4];
		}
	}
}

void Gimli::outputState(uint32 state[]) {
	cout << "State: " << endl;
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 4; c++) {
			cout << hex << setw(9) << state[r + 3 * c] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

uint32 Gimli::getRand32() {
	uint32 sum = rand() % E[16];
	sum = (sum << 16) | (rand() % E[16]);
	sum = sum + previousINT32;
	previousINT32 = sum;
	return sum;
}

uint32 Gimli::getRand16() {
	uint32 sum = getRand32()&0xffff;
	return sum;
}

uint32 Gimli::getRand8() {
	uint32 sum = getRand32() & 0xff;
	return sum;
}

void Gimli::toBit(uint32 s, bool b[], int size) {
	for (int i = 0; i < size; i++) {
		b[i] = (s >> i) & 0x1;
	}
}

uint32 Gimli::toUint32(bool b[], int size) {
	uint32 val = 0;
	for (int i = 0; i < size; i++) {
		if (b[i]) {
			val = val | E[i];
		}
	}
	return val;
}

bool Gimli::BIT(uint32 x, int n) {
	return (x >> n) & 0x1;
}

uint32 Gimli::RL(uint32 word, int n, int wordSize) {
	//n<wordSize
	uint32 mask=0xffffffff;
	if (wordSize < 32) {
		mask = E[wordSize] - 1;
	}
	uint32 high = (word << n) & mask;
	uint32 low = word >> (wordSize - n);
	return high | low;
}

uint32 Gimli::L(uint32 word, int n, int wordSize) {
	uint32 mask = 0xffffffff;
	if (wordSize < 32) {
		mask = E[wordSize] - 1;
	}
	uint32 high = (word << n) & mask;
	return high;
}

void Gimli::SP(uint32 ix, uint32 iy, uint32 iz, uint32 &ox, uint32 &oy, uint32 &oz, int wordSize) {
	int sx = wordSize * 24 / 32;
	int sy = (wordSize * 8 / 32 + 1) % wordSize;
	//ix = ix <<< (s1)
	//iy = iy <<< (s2)
	uint32 tx, ty, tz;
	tx = RL(ix, sx, wordSize);
	ty = RL(iy, sy, wordSize);
	tz = iz;

	//ox = tz ^ ty ^ (tx & ty)<<3
	ox = tz ^ ty ^ L(tx & ty, 3, wordSize);
	//oy = ty ^ tx ^ (tx ^ tz)<<1
	oy = ty ^ tx ^ L(tx | tz, 1, wordSize);
	//oz = tx ^ (tz<<1) ^ (ty & tz)<<2
	oz = tx ^ L(tz, 1, wordSize) ^ L(ty & tz, 2, wordSize);
}

void Gimli::SPInver(uint32 input[],int wordSize) {
	int sx = wordSize * 24 / 32;
	int sy = (wordSize * 8 / 32 + 1) % wordSize;
	uint32 tx = input[0], ty = input[1], tz = input[2];
	bool xBit[32], yBit[32], zBit[32];
	for (int i = 0; i < wordSize; i++) {
		xBit[i] = (tx >> i) & 0x1;
		yBit[i] = (ty >> i) & 0x1;
		zBit[i] = (tz >> i) & 0x1;
	}
	//bit level computation
	bool reXBit[32], reYBit[32], reZBit[32];
	//Firstly, compute reXBit[0],reYBit[0],reZBit[0]
	reXBit[0] = zBit[0];
	reYBit[0] = yBit[0] ^ reXBit[0];
	reZBit[0] = xBit[0] ^ reYBit[0];
	//Secondly, compute reXBit[1],reYBit[1],reZBit[1]
	reXBit[1] = zBit[1] ^ reZBit[0];
	reYBit[1] = yBit[1] ^ reXBit[1] ^ (reXBit[0] | reZBit[0]);
	reZBit[1] = xBit[1] ^ reYBit[1];
	//Thirdly, compute reXBit[2],reYBit[2],reZBit[2]
	reXBit[2] = zBit[2] ^ reZBit[1] ^ (reYBit[0] & reZBit[0]);
	reYBit[2] = yBit[2] ^ reXBit[2] ^ (reXBit[1] | reZBit[1]);
	reZBit[2] = xBit[2] ^ reYBit[2];
	//Recursively recover
	for (int i = 3; i < wordSize; i++) {
		reXBit[i] = zBit[i] ^ reZBit[i - 1] ^ (reYBit[i - 2] & reZBit[i - 2]);
		reYBit[i] = yBit[i] ^ reXBit[i] ^ (reXBit[i - 1] | reZBit[i - 1]);
		reZBit[i] = xBit[i] ^ reYBit[i] ^ (reXBit[i - 3] & reYBit[i - 3]);
	}
	uint32 x = 0;
	uint32 y = 0;
	uint32 z = 0;
	for (int i = 0; i < wordSize; i++) {
		if (reZBit[i]) {
			z |= E[i];
		}
		if (reYBit[i]) {
			y |= E[(i - sy + wordSize) % wordSize];
		}
		if (reXBit[i]) {
			x |= E[(i - sx + wordSize) % wordSize];
		}
	}
	input[0] = x;
	input[1] = y;
	input[2] = z;
}

void Gimli::inverseSPColumn(uint32 state[], int size,int wordSize) {
	uint32 stateTmp[3];
	for (int i = 0; i < size / 3; i++) {
		stateTmp[0] = state[i * 3];
		stateTmp[1] = state[i * 3 + 1];
		stateTmp[2] = state[i * 3 + 2];
		SPInver(stateTmp,wordSize);

		state[i * 3] = stateTmp[0];
		state[i * 3 + 1] = stateTmp[1];
		state[i * 3 + 2] = stateTmp[2];
	}
}

void Gimli::SPColumn(uint32 state[], int size,int wordSize) {
	uint32 x, y, z;
	for (int i = 0; i < size / 3; i++) {
		x = state[3 * i];
		y = state[3 * i + 1];
		z = state[3 * i + 2];
		SP(x, y, z, state[3 * i], state[3 * i + 1], state[3 * i + 2], wordSize);
	}
}

void Gimli::bigSwap(uint32 state[]) {
	swap(state[0], state[6]);
	swap(state[3], state[9]);
}

void Gimli::smallSwap(uint32 state[]) {
	swap(state[0], state[3]);
	swap(state[6], state[9]);
}

void Gimli::check5RPreimageAttack_FindCapacityPart() {
	//To efficiently simulate the 5-round preimage attack,
	//we choose wordSize = 8.
	int wordSize = 8;
	//randomly generate constants
	for (int i = 0; i < 6; i++) {
		constant[i] = getRand8();
	}
	cout << "round constants:" << endl;
	for (int i = 0; i < 6; i++) {
		cout<<hex<<constant[i] << " ";
	}
	cout << endl;
	//Firstly, we choose a challenge
	uint32 hash[2][4];
	for (int i = 0; i < 4; i++) {
		hash[0][i] = rand() % E[8];
		hash[1][i] = rand() % E[8];
	}
	cout << "hash0:" << endl;
	for (int i = 0; i < 4; i++) {
		cout << hex<<hash[0][i] << " ";
	}
	cout << endl;
	cout << "hash1:" << endl;
	for (int i = 0; i < 4; i++) {
		cout << hex << hash[1][i] << " ";
	}
	cout << endl;

	uint32 state[6],s5[12];
	//find the capacity part
	CapHalf *list;
	list = new CapHalf[E[16]];
	//firstly, randomly choose 2^{2wordSize}=2^16 values for s^0[1,2][0,1]
	for (uint32 i = 0; i < E[16]; i++) {
		list[i].in[0] = getRand8();
		list[i].in[1] = getRand8();
		list[i].in[2] = getRand8();
		list[i].in[3] = getRand8();
		state[0] = hash[0][0];
		state[1] = list[i].in[0];
		state[2] = list[i].in[1];
		state[3] = hash[0][1];
		state[4] = list[i].in[2];
		state[5] = list[i].in[3];
		SPColumn(state, 6, wordSize);
		swap(state[0], state[3]);
		state[0] ^= constant[0];
		SPColumn(state, 6, wordSize);
		SPColumn(state, 6, wordSize);
		list[i].row[0] = state[0];
		list[i].row[1] = state[3];

		list[i].row[2] = state[1];
		list[i].row[2] = list[i].row[2]|(state[2]<<8);
		list[i].row[2] = list[i].row[2] | (state[4] << 16);
		list[i].row[2] = list[i].row[2] | (state[5] << 24);
	}
	//sort s
	sort(list, list + E[16], compare);
	//secondly, randomly choose values for s^5[1,2][0,1] until we find a match
	bool find = false;
	uint32 v = 0;
	CapHalf tmpcap;
	uint32 runningTimes = 0;
	while (!find) {
		state[0] = hash[1][0];
		state[1] = getRand8();
		state[2] = getRand8();
		state[3] = hash[1][1];
		state[4] = getRand8();
		state[5] = getRand8();
		s5[0] = hash[1][0];
		s5[1] = state[1];
		s5[2] = state[2];
		s5[3] = hash[1][1];
		s5[4] = state[4];
		s5[5] = state[5];

		state[0] ^= constant[1];
		swap(state[0], state[3]);
		inverseSPColumn(state, 6, wordSize);
		inverseSPColumn(state, 6, wordSize);
		//compare state[1,2,4,5] with list[i].row[2]
		v = state[1];
		v = v | (state[2] << 8);
		v = v | (state[4] << 16);
		v = v | (state[5] << 24);
		tmpcap.row[2] = v;
		find = binary_search(list, list + E[16], tmpcap, compare);
		runningTimes++;
		if (find) {
			int index = binary_search_find_index(list, E[16], tmpcap);
			//cout << list[index].row[0] << " " << list[index].row[1] << endl;
			find=fillLeftWords_FindCap(list,index,state,hash,s5);
		}
		if (find) {
			cout << "running times:" << hex << runningTimes<<" = "<<(double)(runningTimes)/E[16]<<" * 2^16." << endl;
			break;
		}
	}

	delete[]list;
}

bool Gimli::fillLeftWords_FindCap(CapHalf* list, int index, uint32 state[],uint32 hash[][4],uint32 s5filled[]) {
	int wordSize = 8;
	uint32 s3Row[4];
	s3Row[0] = state[0];
	s3Row[1] = state[3];
	s3Row[2] = list[index].row[0];
	s3Row[3] = list[index].row[1];
	uint32 s5[3],sh[6];
	uint32 s0[6];
	CapHalf* capHalf;
	capHalf = new CapHalf[E[7]];
	int cnt = 0;
	for (int i = 0; i < E[16]; i++) {
		s5[0] = hash[1][3];//due to small-swap
		s5[1] = i & 0xff;
		s5[2] = (i >> 8) & 0xff;
		inverseSPColumn(s5, 3, wordSize);
		inverseSPColumn(s5, 3, wordSize);
		//check s5[0]=s3Row[2]
		if (s5[0] == s3Row[2]) {
			//move backwards further
			s5[0] = s3Row[0];//due to big-swap
			inverseSPColumn(s5,3,wordSize);
			inverseSPColumn(s5, 3, wordSize);
			sh[1] = s5[1];
			sh[2] = s5[2];
			sh[3] = s5[0];//due to small-swap
			//check the validity of (sh[1],sh[2],hash[0][2])
			//cout << hex<<hash[0][2] << " " << sh[1] << " " << sh[2] << endl;
			if (property3(hash[0][2], s0[1], s0[2], sh[0], sh[1], sh[2])) {
				//store (s0[1,2], sh[0],s5[0])
				capHalf[cnt].in[0] = i & 0xff;
				capHalf[cnt].in[1] = (i>>8) & 0xff;

				//cout << "sh[0]:" << sh[0]<<" "<<sh[3] << endl;

				capHalf[cnt].row[2] = sh[0] & 0x7f;//only interested in the lower 7 bits
				capHalf[cnt].row[2] = capHalf[cnt].row[2] | (sh[3] << 8);
				cnt++;
				//cout << hex << cnt << endl;
			}
		}
		if (cnt == E[7]) {
			break;
		}
	}
	//sort s
	sort(capHalf, capHalf + cnt, compare);

	//secondly, randomly choose values for s^5[1,2][0,1] until we find a match
	bool find = false;
	uint32 v = 0;
	CapHalf tmpcap;
	for (int i = 0; i < E[16]; i++) {
		s5[0] = hash[1][2];//due to small-swap
		s5[1] = i & 0xff;
		s5[2] = (i >> 8) & 0xff;
		inverseSPColumn(s5, 3, wordSize);
		inverseSPColumn(s5, 3, wordSize);
		//check s5[0]=s3Row[3]
		if (s5[0] == s3Row[3]) {
			//move backwards further
			s5[0] = s3Row[1];//due to big-swap
			inverseSPColumn(s5, 3, wordSize);
			inverseSPColumn(s5, 3, wordSize);
			sh[4] = s5[1];
			sh[5] = s5[2];
			sh[0] = s5[0];//due to small-swap
			//check the validity of (sh[1],sh[2],hash[0][2])
			if (property3(hash[0][3], s0[4], s0[5], sh[3], sh[4], sh[5])) {
				sh[0] = sh[0] & 0x7f;//only interested in the lower 7 bits
				v = sh[0];
				v = v | (sh[3] << 8);
				tmpcap.row[2] = v;
				find = binary_search(capHalf, capHalf + cnt, tmpcap, compare);
				if (find) {
					int index = binary_search_find_index(capHalf, cnt, tmpcap);
					s5filled[6] = hash[1][2];
					s5filled[7] = capHalf[index].in[0];
					s5filled[8] = capHalf[index].in[1];
					s5filled[9] = hash[1][3];
					s5filled[10] = i & 0xff;
					s5filled[11] = (i >> 8) & 0xff;
					cout << "find s5:" << endl;
					outputState(s5filled);

					s5filled[0] ^= constant[1];
					smallSwap(s5filled);
					inverseSPColumn(s5filled,12,wordSize);

					inverseSPColumn(s5filled,12, wordSize);

					bigSwap(s5filled);
					inverseSPColumn(s5filled, 12, wordSize);

					inverseSPColumn(s5filled, 12, wordSize);
					s5filled[0] ^= constant[0];
					smallSwap(s5filled);
					inverseSPColumn(s5filled, 12, wordSize);
					cout << "inverse s5:" << endl;
					outputState(s5filled);
					if (s5filled[0] == hash[0][0] &&
						s5filled[3] == hash[0][1] &&
						s5filled[6] == hash[0][2] &&
						s5filled[9] == hash[0][3]) {
						cout << "match the hash value!" << endl;
					}
					else {
						cout << "not match the hash value!" << endl;
					}
				}
			}
		}
		if (find) {
			break;
		}
	}
	delete[]capHalf;

	if (!find) {
		return false;
	}
	return true;
}

bool Gimli::property3(uint32 ix, uint32& iy, uint32& iz, uint32& ox, uint32 oy, uint32 oz) {
	int wordSize = 8;
	//wordSize = 8
	//ix = ix<<6
	//iy = iy<<3
	
	//oz = ix ^ iz<<1 ^ (iy&iz)<<2
	//oy = iy ^ ix ^ (ix | iz)<<1
	//ox = iz ^ iy ^ (ix&iy)<<3
	if (BIT(ix, 2) != BIT(oz, 0)) {
		//cout << "false" << endl;
		return false;
	}
	//we know ix,oy,oz
	bool IX[8], IY[8], IZ[8], OX[8], OY[8], OZ[8];

	ix = RL(ix, 6, wordSize);
	toBit(ix, IX,8);
	toBit(oy, OY,8);
	toBit(oz, OZ,8);

	IY[0] = OY[0] ^ IX[0];
	IZ[0] = OZ[1] ^ IX[1];

	IY[1] = OY[1] ^ IX[1] ^ (IX[0] | IZ[0]);
	IZ[1] = OZ[2] ^ IX[2] ^ (IY[0] & IZ[0]);
	for (int i = 2; i < 8; i++) {
		IY[i] = OY[i] ^ IX[i] ^ (IX[i - 1] | IZ[i - 1]);
		if(i<7)
			IZ[i] = OZ[i+1] ^ IX[i+1] ^ (IY[i - 1] & IZ[i - 1]);
	}
	//IZ[7]=0
	IZ[7] = 0;//assume
	for (int i = 0; i < 3; i++) {
		OX[i] = IZ[i] ^ IY[i];
	}
	for (int i = 3; i < 8; i++) {
		OX[i] = IZ[i] ^ IY[i] ^ (IX[i - 3] & IY[i - 3]);
	}
	ox=toUint32(OX,8);
	iy = toUint32(IY, 8);
	iz = toUint32(IZ, 8);

	//cout << "true" << endl;
	//cout << hex << "ox: " << ox << endl;
	return true;
}

bool Gimli::check5RPreimageAttack_MatchCapacityPart() {
	int wordSize=8;
	//generate round constants
	for (int i = 0; i < 6; i++) {
		constant[i] = getRand8();
	}
	cout << "round constants:" << endl;
	for (int i = 0; i < 6; i++) {
		cout << hex << constant[i] << " ";
	}
	cout << endl;

	//a random capacity part
	uint32 cap[8];
	for (int i = 0; i < 8; i++) {
		//cap[i] = getRand8();
		cap[i] = 0;//set as all 0
	}
	cout << "challenge capacity part:" << endl;
	cout << hex << cap[0] << " " << cap[1] << " " << cap[2] << " " << cap[3] << endl;
	cout << hex << cap[4] << " " << cap[5] << " " << cap[6] << " " << cap[7] << endl;
	

	//construct tables to store the mapping
	vector<vector<CapHalf> > tableL,tableR;
	tableL.clear();
	tableL.resize(E[16]);
	tableR.clear();
	tableR.resize(E[16]);
	for (int i = 0; i < E[16]; i++) {
		tableL[i].clear();
		tableR[i].clear();
	}
	//construct tables
	constructMappingTable(tableL, constant, 1);
	constructMappingTable(tableR, constant, 0);

	//consider the second message words
	uint32 s5L[6],s5R[6];
	uint32 shL[4], shR[4];//the first row of s^0.5 computed from s5L and s5R, resp.
	CapHalf *list;
	list = new CapHalf[E[17]];//store shL
	int cnt = 0;
	uint32 t[3];
	//construct list
	bool isFull = false;
	for (uint32 i = 0; !isFull && i < E[16]; i++) {
		s5L[0] = i & 0xff;
		s5L[1] = cap[0];
		s5L[2] = cap[4];
		s5L[3] = (i >> 8) & 0xff;
		s5L[4] = cap[2];
		s5L[5] = cap[6];
		//add constants and small swap have been considered
		inverseSPColumn(s5L, 6, wordSize);
		inverseSPColumn(s5L, 6, wordSize);
		swap(s5L[0], s5L[3]);
		inverseSPColumn(s5L, 6, wordSize);
		inverseSPColumn(s5L, 6, wordSize);
		shL[1] = s5L[0]^constant[0];//constant addition
		shL[3] = s5L[3];
		//compute shL[0] and shL[2] using tableL
		//guess shL[0]
		for (uint32 k = 0; !isFull && k < E[8]; k++) {
			t[0] = k;
			t[1] = s5L[1];
			t[2] = s5L[2];
			SPInver(t, wordSize);
			uint32 index = t[1] | (t[2] << 8);
			for (int it = 0; it < !isFull && tableL[index].size(); it++) {
				if (property4(tableL[index][it].row[0], tableL[index][it].row[1], s5L[4], s5L[5], shL[2])) {
					shL[0] = k;//all shL[0,1,2,3] are fixed
					
					list[cnt].in[0] = tableL[index][it].in[0];
					list[cnt].in[1] = tableL[index][it].in[1];
					list[cnt].row[0] = i&0xff;
					list[cnt].row[1] = (i >> 8) & 0xff;

					list[cnt].row[2] = shL[0];
					list[cnt].row[2] = list[cnt].row[2] | (shL[1]<<8);
					list[cnt].row[2] = list[cnt].row[2] | (shL[2] << 16);
					list[cnt].row[2] = list[cnt].row[2] | (shL[3] << 24);
					cnt++;
					if (cnt == E[17]) {
						isFull = true;
						break;
					}
				}
			}
		}
	}
	//search a match
	sort(list, list + cnt, compare);
	//cout << "sorting is over" << endl;
	//cout << "cnt: " << hex<<cnt << endl;

	//match 
	bool isFind = false;
	CapHalf ch;
	uint32 candi = 0;
	for (uint32 i = 0; !isFind && i < E[16]; i++) {
		s5R[0] = i & 0xff;
		s5R[1] = cap[1];
		s5R[2] = cap[5];
		s5R[3] = (i >> 8) & 0xff;
		s5R[4] = cap[3];
		s5R[5] = cap[7];

		//add constants
		s5R[0] ^= constant[1];

		inverseSPColumn(s5R, 6, wordSize);
		inverseSPColumn(s5R, 6, wordSize);
		swap(s5R[0], s5R[3]);
		inverseSPColumn(s5R, 6, wordSize);
		inverseSPColumn(s5R, 6, wordSize);
		shR[0] = s5R[0];
		shR[2] = s5R[3];
		//compute shL[1] and shL[3] using tableR
		//guess shL[0]
		for (uint32 k = 0; !isFind && k < E[8]; k++) {
			t[0] = k;
			t[1] = s5R[1];
			t[2] = s5R[2];
			SPInver(t, wordSize);

			uint32 index = t[1] | (t[2] << 8);
			for (int it = 0; it < !isFull && tableR[index].size(); it++) {
				if (property4(tableR[index][it].row[0], tableR[index][it].row[1], s5R[4], s5R[5], shR[3])) {
					candi++;
					shR[1] = k;//all shR[0,1,2,3] are fixed

					ch.row[2] = shR[0];
					ch.row[2] = ch.row[2] | (shR[1] << 8);
					ch.row[2] = ch.row[2] | (shR[2] << 16);
					ch.row[2] = ch.row[2] | (shR[3] << 24);

					isFind = binary_search(list, list + cnt, ch,compare);
					if (isFind) {
						int in = binary_search_find_index(list, cnt, ch);
						/*//cout << list[in].in[0] << " ";
						cout << tableR[index][it].in[0] << " ";
						cout << list[in].in[1] << " ";
						cout << tableR[index][it].in[1] << " " << endl;

						cout << (i & 0xff) << " ";
						cout << list[in].row[0] << " ";
						cout << ((i >> 8) & 0xff) << " ";
						cout << list[in].row[1] << " " << endl;
						*/

						cout << "find" << endl;

						//verify
						uint32 m[12],n[12];
						for (int i = 0; i < 12; i++) {
							m[i] = 0;
						}
						m[0] = list[in].in[0];
						m[3] = tableR[index][it].in[0];
						m[6] = list[in].in[1];
						m[9] = tableR[index][it].in[1];

						uint32 m0[4];
						m0[0] = m[0];
						m0[1] = m[3];
						m0[2] = m[6];
						m0[3] = m[9];
						cout << "m0: " << m0[0] << " " << m0[1] << " " << m0[2] << " " << m0[3] << endl;

						permutation5R(m, constant, wordSize);
						//outputState(m);

						n[0] = i & 0xff;
						n[1] = cap[0];
						n[2] = cap[4];
						n[3] = list[in].row[0];
						n[4] = cap[1];
						n[5] = cap[5];
						n[6] = (i >> 8) & 0xff;
						n[7] = cap[2];
						n[8] = cap[6];
						n[9] = list[in].row[1];
						n[10] = cap[3];
						n[11] = cap[7];

						inversePermutation5R(n, constant, wordSize);
						//outputState(n);

						uint32 m1[4];
						m1[0] = n[0] ^ m[0];
						m1[1] = n[3] ^ m[3];
						m1[2] = n[6] ^ m[6];
						m1[3] = n[9] ^ m[9];
						cout << "m1: " << m1[0] << " " << m1[1] << " " << m1[2] << " " << m1[3] << endl;
						
						uint32 iv[12] = { 0 };
						iv[0] ^= m0[0];
						iv[3] ^= m0[1];
						iv[6] ^= m0[2];
						iv[9] ^= m0[3];
						permutation5R(iv,constant,wordSize);
						iv[0] ^= m1[0];
						iv[3] ^= m1[1];
						iv[6] ^= m1[2];
						iv[9] ^= m1[3];
						permutation5R(iv, constant, wordSize);
						cout << "after compressing m0 and m1:" << endl;
						outputState(iv);
					}

				}
			}
		}
	}
	//cout << "candi: " << candi << endl;

	for (int i = 0; i < E[16]; i++) {
		tableL[i].clear();
		tableR[i].clear();
	}
	tableL.clear();
	tableR.clear();
	delete[]list;
	return isFind;
}

void Gimli::permutation5R(uint32 s[], uint32 constant[],int wordSize) {
	SPColumn(s, 12, wordSize);
	smallSwap(s);
	s[0] ^= constant[0];
	SPColumn(s, 12, wordSize);
	SPColumn(s, 12, wordSize);
	bigSwap(s);
	SPColumn(s, 12, wordSize);
	SPColumn(s, 12, wordSize);
	smallSwap(s);
	s[0] = s[0] ^ constant[1];
}

void Gimli::inversePermutation5R(uint32 s[], uint32 constant[], int wordSize) {
	s[0] ^= constant[1];
	smallSwap(s);
	inverseSPColumn(s, 12, wordSize);
	inverseSPColumn(s, 12, wordSize);
	bigSwap(s);
	inverseSPColumn(s, 12, wordSize);
	inverseSPColumn(s, 12, wordSize);
	s[0] ^= constant[0];
	smallSwap(s);
	inverseSPColumn(s, 12, wordSize);
}

void Gimli::constructMappingTable(vector<vector<CapHalf> >& table, uint32 constant[], int num) {
	int wordSize = 8;
	//round 1:
	uint32 s[6] = { 0,0,0,0,0,0 };
	for (int i = 0; i < E[16]; i++) {
		s[0] = i & 0xff;
		s[1] = 0;
		s[2] = 0;
		s[3] = (i >> 8) & 0xff;
		s[4] = 0;
		s[5] = 0;
		//SP
		SPColumn(s, 6, wordSize);
		
		//small swap has no influence as the first row are all zero
		if (num == 1) {//there is constant addition
			s[0] ^= constant[0];
		}
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		swap(s[0], s[3]);
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		//store s[1] s[2] s[4] s[5] 
		uint32 index = s[1] | (s[2] << 8);
		CapHalf c;
		c.in[0] = i & 0xff;
		c.in[1] = (i >> 8) & 0xff;
		c.row[0] = s[4];
		c.row[1] = s[5];
		table[index].push_back(c);
	}
}

bool Gimli::property4(uint32 iy, uint32 iz, uint32 oy, uint32 oz, uint32& ox) {
	int wordSize = 8;
	int sx = 6;
	int sy = 3;
	uint32 ix,oyt;
	//iy,iz,oy,oz are known, compute ox

	//x = z ^ y ^ (x & y)<<3
	//y = y ^ x ^ (x | z)<<1
	//z = x ^ (z<<1) ^ (y & z)<<2

	iy = RL(iy, sy, wordSize);
	ix = oz ^ L(iz, 1, wordSize) ^ L(iy & iz, 2, wordSize);
	oyt = iy ^ ix ^ L(ix | iz, 1, wordSize);
	if (oyt != oy) {
		return false;
	}
	ox = iz ^ iy ^ L(ix & iy, 3, wordSize);

	//verify
	uint32 t[3];
	t[0] = ox;
	t[1] = oy;
	t[2] = oz;
	SPInver(t, wordSize);
	t[1] = RL(t[1], sy, wordSize);
	if (t[1] != iy || t[2] != iz) {
		cout << "wrong" << endl;
	}

	return true;
}

bool Gimli::check9RPreimageAttackXOF_FindRate() {
	//wordSize = 8
	int wordSize = 8;
	//generate round constants
	for (int i = 0; i < 6; i++) {
		constant[i] = getRand8();
	}
	cout << "round constants:" << endl;
	for (int i = 0; i < 6; i++) {
		cout << hex << constant[i] << " ";
	}
	cout << endl;
	//hash challenge
	uint32 hash[4];
	for (int i = 0; i < 4; i++) {
		hash[i] = getRand8();
	}

	CapHalf* t0, *t1, *t2, *t3;
	int tSize[4] = { E[16],E[16],E[16],E[16] };
	t0 = new CapHalf[tSize[0]];
	t1 = new CapHalf[tSize[1]];
	t2 = new CapHalf[tSize[2]];
	t3 = new CapHalf[tSize[3]];
	twoRoundsSP(hash[0], t0, tSize[0]);
	twoRoundsSP(hash[1], t1, tSize[1]);
	twoRoundsSP(hash[2], t2, tSize[2]);
	twoRoundsSP(hash[3], t3, tSize[3]);

	uint32 row3[4];
	uint32 row2[4];
	uint32 row1[4];
	uint32 s[6];
	for (int i = 0; i < 4; i++) {
		row2[i] = 0;
		row3[i] = getRand8();
		row1[i] = row3[i];//property 1 
	}
	uint32 r[4];
	CapHalf* listMiddle;
	int size = E[10];
	listMiddle = new CapHalf[size];
	int finalSize = 0;
	CapHalf* listRow;
	int sizeRow = E[16];
	listRow = new CapHalf[sizeRow];
	int cnt = 0;
	bool isMiddleMatched;
	for (int i = 0; i < E[16]; i++) {
		s[0] = i & 0xff;
		s[3] = (i >> 8) & 0xff;
		s[1] = row2[0];
		s[2] = row3[0];
		s[4] = row2[2];
		s[5] = row3[2];
		//SP
		SPColumn(s, 6, wordSize);
		//small swap
		s[0] = row1[1];
		s[3] = row1[3];
		//constant
		s[0] ^= constant[0];
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		swap(s[0], s[3]);
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		r[1] = s[0];
		r[3] = s[3];
		//compute r[0]
		fourRoundsOneColumn(s[1], s[2], hash[0], t0, tSize[0],listMiddle,size,finalSize, false,r[0],r[2]);
		sort(listMiddle, listMiddle + finalSize, compare);
		//compute r[2]
		isMiddleMatched=fourRoundsOneColumn(s[4], s[5], hash[2], t2, tSize[2],listMiddle, finalSize, finalSize, true,r[0],r[2]);
		//test r[0], r[2]
		//cout << hex<<finalSize <<" "<< isMiddleMatched << endl;
		//system("pause");
		if (isMiddleMatched) {
			//cout << "matched in the middle" << endl;
			/*s[0] = r[0];
			s[3] = r[2];
			SPColumn(s, 6, wordSize);
			SPColumn(s, 6, wordSize);
			swap(s[0], s[3]);
			SPColumn(s, 6, wordSize);
			SPColumn(s, 6, wordSize);
			cout << s[0] << " " << s[3] << endl;
			cout << hash[0] << " " << hash[2] << endl;
			//system("pause");*/

			//store r[0],r[1],r[2],r[3]
			listRow[cnt].in[0] = i & 0xff;
			listRow[cnt].in[1] = (i >> 8) & 0xff;
			listRow[cnt].row[2] = r[0];
			listRow[cnt].row[2] = listRow[cnt].row[2] | (r[1] << 8);
			listRow[cnt].row[2] = listRow[cnt].row[2] | (r[2] << 16);
			listRow[cnt].row[2] = listRow[cnt].row[2] | (r[3] << 24);
			cnt++;
			if (cnt == sizeRow) {
				break;
			}
		}
		if (cnt == sizeRow) {
			break;
		}
	}
	//sort listRow
	sort(listRow, listRow + cnt, compare);
	//cout << "sort is over" << endl;
	//cout << "count: " <<cnt<< endl;
	//system("pause");

	//clear listMiddle
	int index = 0;
	bool isMatched = false;
	uint32 matchedNum=0;
	for (int i = 0; i < E[16]; i++) {
		s[0] = i & 0xff;
		s[3] = (i >> 8) & 0xff;
		s[1] = row2[1];
		s[2] = row3[1];
		s[4] = row2[3];
		s[5] = row3[3];
		//SP
		SPColumn(s, 6, wordSize);
		//small swap
		s[0] = row1[0];
		s[3] = row1[2];
		
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		swap(s[0], s[3]);
		SPColumn(s, 6, wordSize);
		SPColumn(s, 6, wordSize);
		r[0] = s[0]^constant[1];//constant addition
		r[2] = s[3];
		//compute r[1]
		fourRoundsOneColumn(s[1], s[2], hash[1], t1, tSize[1], listMiddle, size, finalSize, false, r[1], r[3]);
		sort(listMiddle, listMiddle + finalSize, compare);
		//compute r[3]
		isMiddleMatched = fourRoundsOneColumn(s[4], s[5], hash[3], t3, tSize[3], listMiddle, finalSize, finalSize, true, r[1], r[3]);
		//test r[1], r[3]
		if (isMiddleMatched) {
			//cout << "matched in the middle" << endl;
			/*s[0] = r[1];
			s[3] = r[3];
			SPColumn(s, 6, wordSize);
			SPColumn(s, 6, wordSize);
			swap(s[0], s[3]);
			SPColumn(s, 6, wordSize);
			SPColumn(s, 6, wordSize);
			cout << s[0] << " " << s[3] << endl;
			cout << hash[1] << " " << hash[3] << endl;
			//system("pause");
			*/

			//store r[0],r[1],r[2],r[3]
			CapHalf ch;
			ch.row[2] = r[0];
			ch.row[2] = ch.row[2] | (r[1] << 8);
			ch.row[2] = ch.row[2] | (r[2] << 16);
			ch.row[2] = ch.row[2] | (r[3] << 24);
			//search
			isMatched = binary_search(listRow, listRow + cnt, ch, compare);
			if (isMatched) {
				index = binary_search_find_index(listRow, cnt, ch);
				/*cout << "final found" << endl;
				cout << listRow[index].in[0] << " ";
				cout << (i & 0xff) << " ";
				cout<< listRow[index].in[1] << " ";
				cout << ((i >> 8) & 0xff) << endl;*/

				uint32 res[12];
				res[0] = listRow[index].in[0];
				res[3] = (i & 0xff);
				res[6] = listRow[index].in[1];
				res[9] = ((i >> 8) & 0xff);

				res[1] = 0;
				res[4] = 0;
				res[7] = 0;
				res[10] = 0;

				res[2] = row3[0];
				res[5] = row3[1];
				res[8] = row3[2];
				res[11] = row3[3];

				cout << "found!" << endl;

				cout << "hash challenge: " << endl;
				cout << hash[0] << " " << hash[1] << " " << hash[2] << " " << hash[3] << endl;

				cout << "input:" << endl;
				outputState(res);

				permutation8HalfR(res, constant, wordSize);
				cout << "output after 8.5 rounds:" << endl;
				outputState(res);

				if (res[0] == hash[0] &&
					res[3] == hash[1] &&
					res[6] == hash[2] &&
					res[9] == hash[3]) {
					cout << "match hash challenge!" << endl;
				}
				else {
					cout << "Not match!" << endl;
				}
				return isMatched;
			}
		}
	}

	delete[]listMiddle;
	delete[]t0;
	delete[]t1;
	delete[]t2;
	delete[]t3;
	return isMatched;
}

bool Gimli::fourRoundsOneColumn(uint32 y, uint32 z, uint32 outputX, CapHalf* table, int tableSize, CapHalf* list, int size, int& finalSize, bool isSearch, uint32& r0, uint32& r1) {
	int index = 0;
	bool find = true;
	uint32 s[3],r[2];
	int wordSize = 8;
	CapHalf ch,v;
	int cnt = 0;
	bool isFull = false;
	bool findList = false;
	int listIndex = 0;
	uint32 findNum = 0;
	uint32 searchNum = 0;
	for (int i = 0; !isFull && i < E[8]; i++) {
		s[0] = i;
		s[1] = y;
		s[2] = z;
		SPColumn(s, 3, wordSize);
		SPColumn(s, 3, wordSize);
		ch.row[2] = outputX;
		ch.row[2] = ch.row[2] | (s[1] << 8);
		ch.row[2] = ch.row[2] | (s[2] << 16);
		find = binary_search(table, table + tableSize, ch, compare);

		if (find) {
			findNum++;
			//cout << "find0 is find" << endl;
			index = binary_search_find_index(table, tableSize, ch);
			//store table[index].in[0]
			vector<int> vindexset;
			vindexset.clear();
			int start = index - 10;
			int end = index + 10;
			if (start < 0)
				start = 0;
			if (end >= tableSize)
				end = tableSize;
			for (int st = start; st < end; st++) {
				if (table[st].row[2] == ch.row[2]) {
					vindexset.push_back(st);
				}
			}
			if (!isSearch) {
				for (int st = 0; st < vindexset.size(); st++) {
					list[cnt].in[0] = i;
					list[cnt].row[2] = table[vindexset[st]].in[0];
					list[cnt].row[2] = list[cnt].row[2] | (s[0] << 8);
					cnt++;
					finalSize = cnt;
					if (cnt == size) {
						isFull = true;
						break;
					}
				}
			}
			else {//search table
				for (int st = 0; st < vindexset.size(); st++) {
					v.row[2] = s[0];
					v.row[2] = v.row[2] | (table[vindexset[st]].in[0] << 8);
					//cout << v.row[2] << " " << size << endl;
					findList = binary_search(list, list + size, v, compare);
					searchNum++;
					if (findList) {
						//cout << "find" << endl;
						//system("pause");
						listIndex = binary_search_find_index(list, size, v);
						r0 = list[listIndex].in[0];
						r1 = i;
						return findList;
					}
				}
			}
		}
	}
	return findList;
}

void Gimli::twoRoundsSP(uint32 outputX, CapHalf* table, int &size) {
	uint32 s[3];
	int cnt = 0;
	for (int i = 0; i < E[16]; i++) {
		s[0] = outputX;
		s[1] = i & 0xff;
		s[2] = (i >> 8) & 0xff;
		SPInver(s, 8);
		SPInver(s, 8);
		table[cnt].in[0] = s[0];
		table[cnt].row[2] = outputX;
		table[cnt].row[2] = table[cnt].row[2] | (s[1] << 8);
		table[cnt].row[2] = table[cnt].row[2] | (s[2] << 16);
		cnt++;
		if (cnt == size) {
			break;
		}
	}
	sort(table, table + cnt, compare);
	size = cnt;
}

void Gimli::permutation8HalfR(uint32 s[], uint32 constant[], int wordSize) {
	SPColumn(s, 12, wordSize);
	smallSwap(s);
	s[0] ^= constant[0];
	SPColumn(s, 12, wordSize);
	SPColumn(s, 12, wordSize);
	bigSwap(s);
	SPColumn(s, 12, wordSize);

	SPColumn(s, 12, wordSize);
	smallSwap(s);
	s[0] = s[0] ^ constant[1];
	SPColumn(s, 12, wordSize);
	SPColumn(s, 12, wordSize);
	bigSwap(s);
	SPColumn(s, 12, wordSize);

	SPColumn(s, 12, wordSize);
}

//output the test times to find the desirable (input,output) pair
uint32 Gimli::checkFullRoundDis() {
	//in this case, we choose wordSize=16
	//the complexity of this distinguisher is 2^26
	//the generic complexity is 2^32 
	//we constraint 32 bit conditions rather than 48 bit conditions on the input

	int wordSize = 16;
	int sx = wordSize * 24 / 32;
	int sy = (wordSize * 8 / 32 + 1) % wordSize;
	//firstly, we modify the round constants
	for (int i = 0; i < 6; i++) {
		constant[i] = getRand16();
	}
	cout << "random round constants:" << endl;
	for (int i = 0; i < 6; i++) {
		cout << hex << constant[i] << " ";
	}
	cout << endl;
	//randomly choose a value for s^9
	uint32 s9[6],s13[6];
	uint32 state[12],stateBack[12],stateFor[12],input[12],internalState[25][12];
	bool find[2] = { false,false};
	uint32 runTimes=0;
	uint32 testTimesUpper = E[23];
	uint32 cnt = 0;
	while (!find[1]) {
		while (!find[0]) {//deal with conditions on s^13
			s9[0] = getRand16();
			s9[1] = getRand16();
			s9[2] = getRand16();

			s9[3] = s9[0] ^ constant[2];
			s9[4] = s9[1];
			s9[5] = s9[2];

			//SP -> (SP->B_SW) -> SP -> (SP->S_SW->AC)
			SPColumn(s9, 6, wordSize);
			SPColumn(s9, 6, wordSize);
			swap(s9[0], s9[3]);
			SPColumn(s9, 6, wordSize);
			SPColumn(s9, 6, wordSize);
			if (s9[0] == s9[3]) {
				//record partial positions of s13
				state[1] = s9[1];
				state[2] = s9[2];
				state[3] = s9[0];

				state[7] = s9[4];
				state[8] = s9[5];
				state[9] = s9[3];

				for (uint32 j = 0; j < E[wordSize]; j++) {
					s13[0] = j;
					s13[1] = s9[1];
					s13[2] = s9[2];

					s13[3] = s13[0] ^ constant[3];
					s13[4] = s9[4];
					s13[5] = s9[5];

					//assign to state
					state[0] = s13[0];
					state[6] = s13[3];

					//SP -> (SP->B_SW) -> SP -> (SP->S_SW->AC)
					SPColumn(s13, 6, wordSize);
					SPColumn(s13, 6, wordSize);
					swap(s13[0], s13[3]);
					SPColumn(s13, 6, wordSize);
					SPColumn(s13, 6, wordSize);
					if (s13[0] == s13[3]) {
						find[0] = true;
						cnt = 0;
						break;
					}
				}
			}
		}
		//randomly choose values for state[4,5]=state[10,11]
		state[4] = getRand16();
		state[5] = getRand16();
		state[10] = state[4];
		state[11] = state[5];

		for (int i = 0; i < 12; i++) {
			stateBack[i] = state[i];
			stateFor[i] = state[i];
		}

		//state = s13
		//forwards to s24
		//SP-> (SP->B_SW) -> SP -> (SP->S_SW->AC)
		//SP-> (SP->B_SW) -> SP -> (SP->S_SW->AC)
		//SP-> (SP->B_SW) -> SP
		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[14]);

		SPColumn(stateFor, 12, wordSize);
		bigSwap(stateFor);
		copyState(stateFor, internalState[15]);

		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[16]);

		SPColumn(stateFor, 12, wordSize);
		smallSwap(stateFor);
		stateFor[0] = stateFor[0] ^ constant[4];
		copyState(stateFor, internalState[17]);

		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[18]);

		SPColumn(stateFor, 12, wordSize);
		bigSwap(stateFor);
		copyState(stateFor, internalState[19]);

		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[20]);

		SPColumn(stateFor, 12, wordSize);
		smallSwap(stateFor);
		stateFor[0] = stateFor[0] ^ constant[5];
		copyState(stateFor, internalState[21]);

		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[22]);
		SPColumn(stateFor, 12, wordSize);
		bigSwap(stateFor);
		copyState(stateFor, internalState[23]);
		SPColumn(stateFor, 12, wordSize);
		copyState(stateFor, internalState[24]);

		runTimes++;
		if (stateFor[3] == stateFor[9] && stateFor[4] == stateFor[10] && stateFor[5] == stateFor[11]) {
			//cout << hex << runTimes << endl;
			//backwards
			//(SP->S_SW->AC) -> SP -> (SP->B_SW) -> SP
		    //(SP->S_SW->AC) -> SP -> (SP->B_SW) -> SP
		    //(SP->S_SW->AC) -> SP -> (SP->B_SW) -> SP
			//(SP->S_SW->AC)
			stateBack[0] ^= constant[3];
			smallSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);
			bigSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);

			stateBack[0] ^= constant[2];
			smallSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);
			bigSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);

			stateBack[0] ^= constant[1];
			smallSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);
			bigSwap(stateBack);
			inverseSPColumn(stateBack, 12, wordSize);
			inverseSPColumn(stateBack, 12, wordSize);

			//check stateBack[0][0,1,...,9]=stateBack[6][0,1,...,9]
			uint32 x = stateBack[3] & 0x3ff;
			uint32 y = stateBack[9] & 0x3ff;

			if (x == y) {
				cout << "found!" << endl;
				cout << "running times:" << hex << runTimes << " = " << double(runTimes) / E[26] << " * 2^26" << endl << endl;
				cout << "s0.5[0,1,2]: " << stateBack[3]<<" "<<stateBack[1]<<" "<<stateBack[2] << endl;
				cout << "s0.5[6,7,8]: " << stateBack[9]<< " " << stateBack[7] << " " << stateBack[8] <<endl << endl;
				cout << "s24[Col(1)]: " << stateFor[3] << " " << stateFor[4] << " " << stateFor[5] << endl;
				cout << "s24[Col(3)]: " << stateFor[9] << " " << stateFor[10] << " " << stateFor[11] << endl;
				find[1] = true;

				//output and check
				stateBack[0] ^= constant[0];
				smallSwap(stateBack);
				inverseSPColumn(stateBack, 12, wordSize);

				cout << "input:" << endl;
				outputState(stateBack);
				cout << "output after 24 rounds of permutation:" << endl;
				outputState(stateFor);

				//check the number of identical bits
				cout << "11 identical bits in x" << endl;
				for (int i = -1; i < 10; i++) {
					cout << BIT(stateBack[0], (i + wordSize+1-sx) % wordSize);
				}
				cout << endl;
				for (int i = -1; i < 10; i++) {
					cout << BIT(stateBack[6], (i + wordSize + 1 - sx) % wordSize);
				}
				cout << endl;

				cout << "11 identical bits in y" << endl;
				for (int i = -1; i < 10; i++) {
					cout << BIT(stateBack[1], (i + wordSize + 1 - sy) % wordSize);
				}
				cout << endl;
				for (int i = -1; i < 10; i++) {
					cout << BIT(stateBack[7], (i + wordSize + 1 - sy) % wordSize);
				}
				cout << endl;

				cout << "10 identical bits in z" << endl;
				for (int i = 0; i < 10; i++) {
					cout << BIT(stateBack[2], i);
				}
				cout << endl;
				for (int i = 0; i < 10; i++) {
					cout << BIT(stateBack[8], i);
				}
				cout << endl;

				//check whether stateBack = stateFor (test the correctness of code)
				//permutation(stateBack, constant, wordSize);
				//outputState(stateFor);
				//outputState(stateBack);
			}
		}
		if (!find[1]) {
			cnt++;
			if (cnt == testTimesUpper) {
				find[0] = false;
			}
		}
	}
}

bool Gimli::zeroInternalDiffAttack18R() {
	//First consider the following sequence of operations (9 rounds)
	//(SP)->(SP->B_SW)->(SP)->(SP->S_SW->AC)->
	//(SP)->(SP->B_SW)->(SP)->(SP->S_SW->AC)->
	//(SP)
	int wordSize = 32;
	bool zeroDiffBit[6], zeroDiffBitNew[6];

	uint32 state[12], tmpState[12], stateNew[12], tmpStateNew[12];

	state[0] = getRand32();
	state[1] = getRand32();
	state[2] = getRand32();

	state[6] = state[0] ^ constant[2];
	//state[6] = getRand32();
	state[7] = state[1];
	state[8] = state[2];

	state[3] = getRand32();
	state[4] = getRand32();
	state[5] = getRand32();

	state[9] = state[3];
	state[10] = state[4];
	state[11] = state[5];


	for (int i = 0; i < 12; i++) {
		stateNew[i] = state[i];
	}
	stateNew[0] = state[6];
	stateNew[6] = state[0];

	for (int i = 0; i < 12; i++) {
		tmpState[i] = state[i];
		tmpStateNew[i] = stateNew[i];
	}

	//cout << "A^9:" << endl;
	//outputState(tmpState);
	//cout << "B^9:" << endl;
	//outputState(tmpStateNew);

	//inverse (until S^5)
	tmpState[0] ^= constant[2];
	smallSwap(tmpState);
	inverseSPColumn(tmpState, 12,wordSize);

	inverseSPColumn(tmpState, 12,wordSize);

	bigSwap(tmpState);
	inverseSPColumn(tmpState, 12,wordSize);

	inverseSPColumn(tmpState, 12,wordSize);

	//new state
	tmpStateNew[0] ^= constant[2];
	smallSwap(tmpStateNew);
	inverseSPColumn(tmpStateNew, 12,wordSize);

	inverseSPColumn(tmpStateNew, 12,wordSize);

	bigSwap(tmpStateNew);
	inverseSPColumn(tmpStateNew, 12,wordSize);

	inverseSPColumn(tmpStateNew, 12,wordSize);

	//further inverse until S^0
	tmpState[0] ^= constant[1];
	smallSwap(tmpState);
	inverseSPColumn(tmpState, 12,wordSize);

	inverseSPColumn(tmpState, 12,wordSize);

	bigSwap(tmpState);
	inverseSPColumn(tmpState, 12,wordSize);

	inverseSPColumn(tmpState, 12,wordSize);

	tmpState[0] ^= constant[0];
	smallSwap(tmpState);

	//cout << "A^0.5:" << endl;
	//outputState(tmpState);

	inverseSPColumn(tmpState, 12,wordSize);
	//cout << "A^0:" << endl;
	//outputState(tmpState);

	//new state
	tmpStateNew[0] ^= constant[1];
	smallSwap(tmpStateNew);
	inverseSPColumn(tmpStateNew, 12,wordSize);

	inverseSPColumn(tmpStateNew, 12,wordSize);

	bigSwap(tmpStateNew);
	inverseSPColumn(tmpStateNew, 12,wordSize);

	inverseSPColumn(tmpStateNew, 12,wordSize);

	tmpStateNew[0] ^= constant[0];
	smallSwap(tmpStateNew);
	//cout << "B^0.5:" << endl;
	//outputState(tmpStateNew);
	inverseSPColumn(tmpStateNew, 12,wordSize);
	//cout << "B^0:" << endl;
	//outputState(tmpStateNew);

	//distinguisher: tmpState[0][8]=tmpState[6][8]
	//distinguisher: tmpState[1][23]=tmpState[7][23]
	zeroDiffBit[4] = BIT(tmpState[0], 8) ^ BIT(tmpState[6], 8);
	zeroDiffBit[5] = BIT(tmpState[1], 23) ^ BIT(tmpState[7], 23);

	zeroDiffBitNew[4] = BIT(tmpStateNew[0], 8) ^ BIT(tmpStateNew[6], 8);
	zeroDiffBitNew[5] = BIT(tmpStateNew[1], 23) ^ BIT(tmpStateNew[7], 23);

	//forward direction
	for (int i = 0; i < 2; i++) {
		SPColumn(state, 12,wordSize);

		SPColumn(state, 12,wordSize);
		bigSwap(state);

		SPColumn(state, 12,wordSize);

		SPColumn(state, 12,wordSize);
		smallSwap(state);

		state[0] ^= constant[i + 3];
	}
	//cout << "A^17:" << endl;
	//outputState(state);
	SPColumn(state, 12,wordSize);
	//cout << "A^18:" << endl;
	//outputState(state);

	for (int i = 0; i < 2; i++) {
		SPColumn(stateNew, 12,wordSize);

		SPColumn(stateNew, 12,wordSize);
		bigSwap(stateNew);

		SPColumn(stateNew, 12,wordSize);

		SPColumn(stateNew, 12,wordSize);
		smallSwap(stateNew);

		stateNew[0] ^= constant[i + 3];
	}
	//cout << "B^17:" << endl;
	//outputState(stateNew);
	SPColumn(stateNew, 12,wordSize);
	//cout << "B^18:" << endl;
	//outputState(stateNew);

	zeroDiffBit[0] = BIT(state[3], 0) ^ BIT(stateNew[9], 0);
	zeroDiffBit[1] = BIT(state[3], 1) ^ BIT(stateNew[9], 1);
	zeroDiffBit[2] = BIT(state[3], 2) ^ BIT(stateNew[9], 2);
	zeroDiffBit[3] = BIT(state[4], 0) ^ BIT(state[5], 0) ^ BIT(stateNew[10], 0) ^ BIT(stateNew[11], 0);

	zeroDiffBitNew[0] = BIT(stateNew[3], 0) ^ BIT(state[9], 0);
	zeroDiffBitNew[1] = BIT(stateNew[3], 1) ^ BIT(state[9], 1);
	zeroDiffBitNew[2] = BIT(stateNew[3], 2) ^ BIT(state[9], 2);
	zeroDiffBitNew[3] = BIT(stateNew[4], 0) ^ BIT(stateNew[5], 0) ^ BIT(state[10], 0) ^ BIT(state[11], 0);

	for (int i = 0; i < 6; i++) {
		if (zeroDiffBit[i] || zeroDiffBitNew[i]) {
			return false;
		}
	}

	return true;
}



