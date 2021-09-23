#include<iostream>
#include<random>
#include<chrono>
#include<pthread.h>
#include<windows.h>
#include<algorithm>
#include<vector>

using namespace std;

random_device device;
mt19937 gen(device());

class mutation_model {
public:
    long long int num_bases;
    float mutation_probability;
    int k;
    vector<char> mutated;
    bernoulli_distribution bernoulli;

    mutation_model(long long int num_total_bases, float mut_rate, int kmer_size) {
        num_bases = num_total_bases;
        mutated.resize(num_bases);
        k = kmer_size;
        mutation_probability = mut_rate;
        bernoulli_distribution temp(mut_rate);
        bernoulli = temp;
    }

    void generate_bases_trivial() {
        for (long long int i = 0; i < num_bases; i++) {
            bernoulli(gen) ? mutated[i] = 1 : mutated[i] = 0;
        }
    }

    void print_bases() {
        for (long long int i = 0; i < num_bases; i++) {
            cout << (int)mutated[i];
        }
        cout << endl;
    }

    int get_n_mut_trivial() {
        long long int n_mut = 0;
        for (long long int i = 0; i < num_bases - k; i++) {
            for (long long int j = i; j < i+k; j++) {
                if (mutated[j]) {
                    n_mut ++;
                    break;
                }
            }
        }
        return n_mut;
    }

    int get_n_mut_eff() {
        long long int n_mut = 0;
        int sum = 0;
        for (long long int i = 0; i < k; i++) {
            sum += mutated[i];
        }
        if (sum > 0) {
            n_mut ++;
        }
        for (long long int i = 1; i < num_bases - k; i++) {
            sum = sum - mutated[i-1] + mutated[i+k-1];
            if (sum > 0) n_mut++;
        }
        return n_mut;
    }

    void shuffle_bases() {
        shuffle (mutated.begin(), mutated.end(), std::default_random_engine(5635635));
    }

    float get_scaled_jaccard_index(float scale_factor) {
        long long int scaled_length = num_bases * scale_factor;
        long long int num_mutated_kmers = 0;
        int num_total_kmers = 0;
        int sum = 0;
        for (long long int i = 0; i < k; i++) {
            sum += mutated[i];
        }
        if (sum > 0) {
            num_mutated_kmers = 1;
        }
        num_total_kmers = 0;
        for (long long int i = 2; i < scaled_length; i++) {
            sum = sum - mutated[i-1] + mutated[i+k-1];
            if (sum > 0) num_mutated_kmers++;
            num_total_kmers++;
        }
        return 1.0 * (num_total_kmers - num_mutated_kmers) / (num_total_kmers + num_mutated_kmers);
    }

};

int main()
{
    long long int L = 1000000000;
    int k = 21;
    float mut_rate = 0.01;
    float scale_factor = 0.1;

    mutation_model m(L+k-1, mut_rate, k);
    m.generate_bases_trivial();
    m.shuffle_bases();
    cout << m.get_scaled_jaccard_index(scale_factor) << endl;
    //m.shuffle_bases();

    return 0;
}
