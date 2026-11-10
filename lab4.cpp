
	#include <iostream>
	#include <stdlib.h>
	#include <chrono>
	#include<vector>

	const int FIELD_SIZE = 173;
	const int ARRAY_SIZE = FIELD_SIZE / 64 + 1;
	const int POLYNOMIAL_DEGREES[5] = { 0, 1, 2, 10, 173 };



	class ElementF {
	public:
		uint64_t coefficients[FIELD_SIZE / 64 + 1] = { 0 };

		ElementF();

		ElementF(const ElementF& other);

		ElementF(std::string binary);

		ElementF operator + (const ElementF& other);

		ElementF operator * (const ElementF& other);

		ElementF RightCycleShift();

		ElementF LeftCycleShift();

		ElementF square();

		ElementF trace();

		bool operator == (const ElementF& other);

		ElementF power(ElementF& other);

		ElementF inv();

		ElementF operator & (const ElementF& other);

		void set_bit(int index);
		const int get_bit(int index);

		int weight();

		friend std::ostream& operator<<(std::ostream& os, const ElementF& element);
	};

	ElementF* calculateMultiplicativeMatrix();


	ElementF& one();

	/*ElementF& zero();*/

	static ElementF* multiplicativeMatrix = calculateMultiplicativeMatrix();

	ElementF ConvertToNumber(std::string binary);

	std::string ConvertToBinary(const ElementF& element);



	ElementF::ElementF()
	{
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			coefficients[i] = 0;
		}
	}

	ElementF::ElementF(const ElementF& other)
	{
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			coefficients[i] = other.coefficients[i];
		}
	}

	ElementF::ElementF(std::string binary)
	{
		ElementF other = ConvertToNumber(binary);
		*this = ElementF(other);
	}

	ElementF ElementF::operator+(const ElementF& other)
	{
		ElementF result;
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			result.coefficients[i] = coefficients[i] ^ other.coefficients[i];
		}

		return result;
	}

	ElementF ElementF::operator*(const ElementF& other)
	{
		ElementF result;
		ElementF left(*this);
		ElementF right(other);


		for (int i = 0; i < FIELD_SIZE; i++) {
			ElementF temp;
			for (int j = 0; j < FIELD_SIZE; j++) {
				int bit = (left & multiplicativeMatrix[j]).weight() & 1;
				if (bit == 1) {
					temp.set_bit(j);
				}
			}
			temp = temp & right;
			if ((temp.weight() & 1) == 1) {
				result.set_bit(i);
			}
			left = left.LeftCycleShift();
			right = right.LeftCycleShift();
		}

		return result;
	}

	ElementF ElementF::RightCycleShift()
	{
		ElementF result;
		int right = FIELD_SIZE % 64;


		uint64_t temp = coefficients[ARRAY_SIZE - 1] >> (right - 1);
		temp = temp & 1;
		result.coefficients[0] = (coefficients[0] << 1) + temp;
		for (int i = 1; i < ARRAY_SIZE; i++) {
			result.coefficients[i] = (coefficients[i] << 1) + (coefficients[i - 1] >> 63);
		}
		result.coefficients[ARRAY_SIZE - 1] = result.coefficients[ARRAY_SIZE - 1] - (static_cast<uint64_t>(1) << right);
		return result;
	}

	ElementF ElementF::LeftCycleShift()
	{
		ElementF result;
		int right = FIELD_SIZE % 64;


		for (int i = 0; i < ARRAY_SIZE - 1; i++) {
			result.coefficients[i] = (coefficients[i] >> 1) + (coefficients[i + 1] << 63);
		}
		int temp = coefficients[0] & 1;

		result.coefficients[ARRAY_SIZE - 1] = (coefficients[ARRAY_SIZE - 1] >> 1) + (static_cast<uint64_t>(temp) << right - 1);

		return result;
	}

	ElementF ElementF::square()
	{
		return this->RightCycleShift();
	}

	ElementF ElementF::trace()
	{
		if ((weight() & 1) == 1) {
			return ElementF(one());
		}
		return ElementF();

	}

	bool ElementF::operator==(const ElementF& other)
	{
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			if (coefficients[i] == other.coefficients[i]) {
				continue;
			}
			else {
				return false;
			}
		}
		return true;
	}



	ElementF ElementF::operator&(const ElementF& other)
	{
		ElementF result;
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			result.coefficients[i] = coefficients[i] & other.coefficients[i];
		}

		return result;
	}


	void ElementF::set_bit(int index)
	{
		coefficients[index / 64] ^= static_cast<unsigned long long>(1) << (index % 64);
	}

	const int ElementF::get_bit(int index)
	{
		return (coefficients[index / 64] >> (index % 64)) & 1;
	}

	int ElementF::weight()
	{
		int weight = 0;
		for (int i = 0; i <= FIELD_SIZE / 64; i++) {
			weight += _mm_popcnt_u64(coefficients[i]);
		}
		return weight;
	}

	std::ostream& operator<<(std::ostream& os, const ElementF& element)
	{
		std::string res = ConvertToBinary(element);
		os << res;
		return os;
	}

	ElementF* calculateMultiplicativeMatrix()
	{
		static ElementF multmatrix[FIELD_SIZE];
		int p = FIELD_SIZE * 2 + 1;
		int powers[FIELD_SIZE] = { 0 };
		powers[0] = 1;
		bool temp[FIELD_SIZE][FIELD_SIZE] = { false };
		for (int i = 1; i < FIELD_SIZE; i++) {
			powers[i] = 2 * powers[i - 1] % p;
		}
		for (int i = 0; i < FIELD_SIZE; i++) {
			for (int j = i; j < FIELD_SIZE; j++) {

				if ((powers[i] + powers[j]) % p == 1 ||
					(powers[i] + powers[j]) % p == p - 1 ||
					(powers[i] - powers[j]) % p == 1 ||
					(powers[i] - powers[j]) % p == -1) {
					temp[i][j] = true;
					temp[j][i] = true;
				}
			}
		}
		for (int i = 0; i < FIELD_SIZE; i++) {
			ElementF column;
			for (int j = 0; j < FIELD_SIZE; j++) {
				if (temp[i][j]) {
					column.set_bit(j);
				}
			}
			multmatrix[i] = column;
		}
		return multmatrix;
	}


	ElementF ElementF::power(ElementF& exp) {
		ElementF result;
		result.coefficients[0] =( unsigned long long)  (- 1);
		result.coefficients[1] = (unsigned long long)  (-1);
		result.coefficients[2] = ((unsigned long long)  (1) << 45) - 1;

		ElementF base = *this;

		for (int i = FIELD_SIZE - 1; i >= 0; --i) {
			result = result.square();
			if (exp.get_bit(i)) {
				result = result * base;
			}
		}
		return result;
	}



	ElementF ElementF::inv() {
		ElementF temp(*this);
		ElementF result = one();
		for (int i = 0; i < FIELD_SIZE - 1; i++) {
			temp = temp.square();
			result = result * temp;
		}
		return result;
	}








	ElementF& one() {
		ElementF result;
		result.set_bit(0);
		return result;
	}


	




	ElementF ConvertToNumber(std::string binary)
	{
		ElementF result;
		for (int i = 0; i < FIELD_SIZE; i++) {
			if (binary.at(i) == '1') {
				result.coefficients[i / 64] += static_cast<unsigned long long>(1) << (i % 64);
			}
		}
		return result;
	}

	std::string ConvertToBinary(const ElementF& element) {
		std::string binary = "";
		for (int i = 0; i < FIELD_SIZE; i++) {
			int bit = (element.coefficients[i / 64] >> (i % 64)) & 1;
			if (bit == 1) {
				binary += '1';
			}
			else {
				binary += '0';
			}
		}
		return binary;
	}



	ElementF generateElementF() {
		std::string binary = "01";
		std::string result = "";
		for (int i = 0; i < FIELD_SIZE; i++) {
			result += binary[rand() % 2];
		}
		return ElementF(result);
	}

	void measureTime(int numberOfMeasurements) {
		long long nanoseconds = 0;
		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();
			ElementF second = generateElementF();

			const auto start{ std::chrono::high_resolution_clock::now() };
			first + second;
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for addition: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;

		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();
			ElementF second = generateElementF();

			const auto start{ std::chrono::high_resolution_clock::now() };
			first* second;
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for multiplication: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;

		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();


			const auto start{ std::chrono::high_resolution_clock::now() };
			first.square();
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for square: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;

		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();
			ElementF second = generateElementF();

			const auto start{ std::chrono::high_resolution_clock::now() };
			first.power(second);
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for power: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;

		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();
			ElementF second = generateElementF();

			const auto start{ std::chrono::high_resolution_clock::now() };
			first.inv();
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for inverse: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;



		for (int i = 0; i < numberOfMeasurements; i++) {
			ElementF first = generateElementF();
			ElementF second = generateElementF();

			const auto start{ std::chrono::high_resolution_clock::now() };
			first.trace();
			const auto end{ std::chrono::high_resolution_clock::now() };
			const std::chrono::duration<int64_t, std::nano> elapsed_seconds{ end - start };
			nanoseconds += elapsed_seconds.count();
		}
		std::cout << "Elapsed nanoseconds for trace: " << nanoseconds / numberOfMeasurements << '\n';
		nanoseconds = 0;
	}



	int main()
	{
		ElementF f("11000000011010110111100100001101000111001010111111001110000101111010010110110001101000010011100101100111111101010001111001010001101011100101100010100110110000011100000010101");
		ElementF g("10110111101110011000111000000110101101011111110101111011000110110100001011001101101110011010010001011010001001101101100111111100111101111111111011100110100000101010011000000");

		ElementF h("01101010001111010110111001101111001000111110010010011001000101000101000010101000000001100101001011000101011100010011011111001011111110001101011010000101010010111011111100111");


		std::cout << "A = " << f << "\n\n";
		std::cout << "B = " << g << "\n\n";
		std::cout << "N = " << h << "\n\n";


		ElementF add = f + g;
		std::cout << "A+B = " << add << "\n\n";

		ElementF mult = f * g;
		std::cout << "A*B = " << mult << "\n\n";

		ElementF square = f.square();
		std::cout << "A^2 = " << square << "\n\n";

		ElementF inv = f.inv();
		std::cout << "A^-1 = " << inv << "\n\n";

		ElementF pow = f.power(h);
		std::cout << "A^N = " << pow << "\n\n";


		ElementF trace = f.trace();

		std::cout << "Tr(A) = " << trace << "\n\n";

		measureTime(100);
	}


