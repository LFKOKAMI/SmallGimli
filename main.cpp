#include "Gimli.h"
#include <ctime>
#include <iostream>
using namespace std;

void testSP() {
	srand(time(NULL));
	int testTimes;
	int wordSize=8;
	Gimli gimli;
	uint32 input[3], mediate[3];
	for (uint32 i = 0; i < E[24]; i++) {
		input[0] = i & 0xff;
		input[1] = (i >> 8) & 0xff;
		input[2] = (i >> 12) & 0xff;

		//cout << hex << input[0] << " " << input[1] << " " << input[2] << endl;
		gimli.SP(input[0], input[1], input[2], mediate[0], mediate[1], mediate[2], wordSize);
		//cout << hex<<mediate[0]<<" " << mediate[1] <<" "<< mediate[2] << endl;
		gimli.SPInver(mediate, wordSize);

		if (mediate[0] != input[0] || mediate[1] != input[1] || mediate[2] != input[2]) {
			cout << "wrong" << endl;
		}
		//system("pause");
	}
}

void test18RoundDistinguisher() {
	Gimli gimli;
	int testTimes = 10000;
	bool flag = true;
	int successfulTimes = 0;
	cout << "test times:" << testTimes << endl;
	while (testTimes--) {
		if (gimli.zeroInternalDiffAttack18R()) {
			successfulTimes++;
		}
	}
	cout << "successful times:" << successfulTimes << endl;
}

void testFullRoundDistinguisher() {
	Gimli gimli;
	gimli.checkFullRoundDis();
}

void test5RPreimageFindCap() {
	Gimli gimli;
	gimli.check5RPreimageAttack_FindCapacityPart();
}

void test5RPreimageMatchCap() {
	Gimli gimli;
	bool isFind = false;
	int cnt = 0;
	while (!isFind) {
		cnt++;
		cout <<endl<< "times: " << cnt << endl;
		isFind=gimli.check5RPreimageAttack_MatchCapacityPart();
	}
	cout << "total running times: " << cnt << " * 2^24" << endl;
}

void test9RXOF() {
	Gimli gimli;
	bool isFind = false;
	int cnt = 0;
	while (!isFind) {
		cnt++;
		cout << endl << "times: " << cnt << endl;
		isFind = gimli.check9RPreimageAttackXOF_FindRate();
	}
	cout << "total running times: " << cnt << " * 2^24" << endl;
}

int main() {
	cout << "#########################################################" << endl << endl;
	cout << "1 -> test the 18-round distinguisher" << endl << endl;
	cout << "2 -> test the full-round distinguisher" << endl << endl;
	cout << "3 -> test the 5-round preimage attack: find capacity" << endl << endl;
	cout << "4 -> test the 5-round preimage attack: match capacity" << endl << endl;
	cout << "5 -> test the 9-round preimage attack: match rate" << endl << endl;
	cout << "#########################################################" << endl;
	int cmd;
	cout << "please input command:";
	while (cin >> cmd) {
		if (cmd == 1) {
			test18RoundDistinguisher();
		}
		else if (cmd == 2) {
			testFullRoundDistinguisher();
		}
		else if (cmd == 3) {
			test5RPreimageFindCap();
		}
		else if (cmd == 4) {
			test5RPreimageMatchCap();
		}
		else if (cmd == 5) {
			test9RXOF();
		}
		else {
			cout << "wrong command!" << endl;
			break;
		}

		cout << endl;
		cout << "#########################################################" << endl << endl;
		cout << "1 -> test the 18-round distinguisher" << endl << endl;
		cout << "2 -> test the full-round distinguisher" << endl << endl;
		cout << "3 -> test the 5-round preimage attack: find capacity" << endl << endl;
		cout << "4 -> test the 5-round preimage attack: match capacity" << endl << endl;
		cout << "5 -> test the 9-round preimage attack: match rate" << endl << endl;
		cout << "#########################################################" << endl;
		cout << endl<<"please input command:";
	}
	//test18RoundDistinguisher();
	//testFullRoundDistinguisher();
	//test5RPreimageFindCap();
	//test5RPreimageMatchCap();
	//test9RXOF();
	return 1;
}