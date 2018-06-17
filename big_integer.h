#ifndef BIG_INTEGER_H
#define BIG_INTEGER_H

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <algorithm>

struct big_integer {
    big_integer();

    big_integer(big_integer const &other);

    big_integer(int a);

    explicit big_integer(std::string const &str);

    ~big_integer();

    big_integer &operator=(big_integer const &other);

    big_integer &operator+=(big_integer const &rhs);

    big_integer &operator-=(big_integer const &rhs);

    big_integer &operator*=(big_integer const &rhs);

    big_integer &operator/=(big_integer const &rhs);

    big_integer &operator%=(big_integer const &rhs);

    big_integer &operator&=(big_integer const &rhs);

    big_integer &operator|=(big_integer const &rhs);

    big_integer &operator^=(big_integer const &rhs);

    big_integer &operator<<=(int rhs);

    big_integer &operator>>=(int rhs);

    big_integer operator+() const;

    big_integer operator-() const;

    big_integer operator~() const;

    big_integer &operator++();

    big_integer operator++(int);

    big_integer &operator--();

    big_integer operator--(int);

    friend bool operator==(big_integer const &a, big_integer const &b);

    friend bool operator!=(big_integer const &a, big_integer const &b);

    friend bool operator<(big_integer const &a, big_integer const &b);

    friend bool operator>(big_integer const &a, big_integer const &b);

    friend bool operator<=(big_integer const &a, big_integer const &b);

    friend bool operator>=(big_integer const &a, big_integer const &b);

    friend int8_t comp(big_integer const &a, big_integer const &b);

    friend int8_t comp_abs(big_integer const &a, big_integer const &b);

    friend std::string to_string(big_integer const &a);

    void mul_short(const big_integer &a, const unsigned int &x);

private:
    std::vector<unsigned int> data;
    bool sign;

    bool is_zero() const;

    void shrink_zeroes();

    size_t length() const;

    unsigned int get_digit(const size_t &i) const;

    void set_digit(const size_t &i, const unsigned int &val);

    void resize(const size_t &ns, const unsigned int &nv);

    unsigned int back() const;

    void len_prepare(const size_t &new_size);

    void reverse_digits();

    void complem_two();

    void recomplem_two();

    void ar_gen(const big_integer &r, const bool &rev);

    void add_abs(const big_integer &r, const size_t &offset);

    void sub_abs(const big_integer &l, const big_integer &r);

    unsigned int div_short(const unsigned int &a);

    template<class FunctorT>
    big_integer &bit_gen(const big_integer &rhs, FunctorT functor);

    void bit_prepare(const big_integer &r);

    big_integer divide(const big_integer &rhs, const bool &mode);

};

big_integer operator+(big_integer a, big_integer const &b);

big_integer operator-(big_integer a, big_integer const &b);

big_integer operator*(big_integer a, big_integer const &b);

big_integer operator/(big_integer a, big_integer const &b);

big_integer operator%(big_integer a, big_integer const &b);

big_integer operator&(big_integer a, big_integer const &b);

big_integer operator|(big_integer a, big_integer const &b);

big_integer operator^(big_integer a, big_integer const &b);

big_integer operator<<(big_integer a, int b);

big_integer operator>>(big_integer a, int b);

bool operator==(big_integer const &a, big_integer const &b);

bool operator!=(big_integer const &a, big_integer const &b);

bool operator<(big_integer const &a, big_integer const &b);

bool operator>(big_integer const &a, big_integer const &b);

bool operator<=(big_integer const &a, big_integer const &b);

bool operator>=(big_integer const &a, big_integer const &b);

int8_t comp(big_integer const &a, big_integer const &b);

int8_t comp_abs(big_integer const &a, big_integer const &b);

std::string to_string(big_integer const &a);

std::ostream &operator<<(std::ostream &s, big_integer const &a);

#endif // BIG_INTEGER_H
