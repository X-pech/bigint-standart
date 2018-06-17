#include "big_integer.h"

#include <cstring>
#include <stdexcept>

typedef unsigned int ui;
typedef unsigned long long ui64;
typedef long long i64;

const ui MAX_DIGIT = UINT32_MAX;
const ui64 BASE = (1LL << 32);
big_integer ZERO = big_integer(0);
big_integer ONE = big_integer(1);
const ui POWER = 32;
const ui BASE_10 = 1000000000;

template<typename T>
ui64 sui64c(T x) {
    return static_cast<ui64>(x);
}

template<typename T>
ui suic(T x) {
    return static_cast<ui>(x);
}

template<typename T>
i64 si64c(T x) {
    return static_cast<i64>(x);
}

const ui HALFBASE = suic(sui64c(1) << (POWER - 1));

size_t big_integer::length() const {
    return data.size();
}

void big_integer::resize(const size_t &ns, const unsigned int &nv = 0) {
    data.resize(ns, nv);
}

void big_integer::len_prepare(const size_t &new_size) {
    if (length() < new_size) {
        resize(new_size, 0);
    }
}

bool big_integer::is_zero() const {
    return (length() == 1 && get_digit(0) == 0);
}

void big_integer::shrink_zeroes() {
    while (length() && !data.back()) {
        data.pop_back();
    }
    if (!length()) {
        data.push_back(0);
        sign = false;
    }
}

ui big_integer::get_digit(const size_t &i) const {
    return (i < length() ? data[i] : 0);
}

ui big_integer::back() const {
    return data.back();
}

void big_integer::set_digit(const size_t &i, const unsigned int &val) {
    if (i < length()) {
        data[i] = val;
    }
}

void big_integer::reverse_digits() {
    for (ui &i : data) {
        i = ~i;
    }
}

void big_integer::complem_two() {
    bool s = sign;
    sign = false;
    reverse_digits();
    *this += 1;
    sign = s;
}

void big_integer::recomplem_two() {
    bool s = sign;
    sign = false;
    *this -= 1;
    reverse_digits();
    sign = s;
}

int8_t comp(big_integer const &a, big_integer const &b) {
    if (a.is_zero() && b.is_zero()) {
        return 0;
    }

    if (a.sign != b.sign) {
        return (a.sign ? -1 : 1);
    }

    return comp_abs(a, b);
}

int8_t comp_abs(big_integer const &a, big_integer const &b) {
    if (a.length() != b.length()) {
        return (a.length() < b.length() ? -1 : 1);
    }

    size_t m = std::max(a.length(), b.length());
    for (size_t i = 0; i < m; i++) {
        if (a.get_digit(m - i - 1) != b.get_digit(m - i - 1)) {
            return (a.get_digit(m - i - 1) < b.get_digit(m - i - 1) ? -1 : 1);
        }
    }
    return 0;
}

void big_integer::add_abs(const big_integer &r, const size_t &offset) {
    len_prepare(r.length() + offset);
    ui64 cf = 0;
    for (size_t i = 0; i < length(); i++) {
        cf += data[i];
        cf += (i - offset < 0 ? 0 : r.get_digit(i - offset));
        data[i] = suic(cf);
        cf >>= POWER;
    }
    if (cf) {
        data.push_back(suic(cf));
    }
    shrink_zeroes();
}

void big_integer::sub_abs(const big_integer &l, const big_integer &r) {
    i64 s = 0;
    bool cf = false;
    for (size_t i = 0; i < l.length(); i++) {
        s = si64c(l.get_digit(i)) - r.get_digit(i) - si64c(cf);
        if (s < 0) {
            cf = true;
            s += BASE;
        } else {
            cf = false;
        }
        data[i] = suic(s);
    }
}

void big_integer::ar_gen(const big_integer &r, const bool &rev) {
    if (r == ZERO) {
        return;
    }

    if ((sign == r.sign) ^ rev) {
        add_abs(r, 0);
        shrink_zeroes();
        return;
    }

    int8_t cr = comp_abs(*this, r);

    if (cr == 0) {
        *this = ZERO;
        return;
    }

    if (cr == -1) {
        sub_abs(r, *this);
        sign = r.sign ^ rev;
    } else {
        sub_abs(*this, r);
    }
    shrink_zeroes();
}

big_integer::big_integer()
        : data({0}), sign(false) {}

big_integer::big_integer(big_integer const &other) = default;

big_integer::big_integer(int a)
        : data({suic((i64) a * ((a < 0) ? -1LL : 1LL))}), sign(a < 0) {}

big_integer::big_integer(std::string const &str)
        : data({0}), sign(false) {
    bool rs = (str[0] == '-');

    ui cd = 0, local_base_10 = 1;
    for (auto i = str.begin() + (rs); i != str.end(); i++) {
        cd = 0;
        local_base_10 = 1;
        for (size_t j = 0; j < 9 && i != str.end(); i++, j++) {
            if (*i < '0' || '9' < *i) {
                throw std::runtime_error("Invalid char!");
            }
            cd *= 10;
            cd += suic(*i - '0');
            local_base_10 *= 10;
        }
        i--;
        mul_short(*this, local_base_10);
        *this += cd;
    }
    sign = rs;
    shrink_zeroes();
}


big_integer::~big_integer() {
    data = std::vector<ui>();
    sign = false;
}

big_integer &big_integer::operator=(big_integer const &other)= default;

big_integer &big_integer::operator+=(big_integer const &rhs) {
    ar_gen(rhs, false);
    return *this;
}

big_integer &big_integer::operator-=(big_integer const &rhs) {
    ar_gen(rhs, true);
    return *this;
}

void big_integer::mul_short(const big_integer &a, const unsigned int &x) {
    size_t alo = a.length();
    len_prepare(alo + 3);

    ui64 cf = 0;
    for (size_t i = 0; i < alo; i++) {
        cf = sui64c(a.get_digit(i)) * sui64c(x) + cf;
        data[i] = suic(cf);
        cf >>= POWER;
    }
    if (cf) {
        data[alo] = suic(cf);
    } else {
        data[alo] = 0;
    }
    shrink_zeroes();
}

big_integer &big_integer::operator*=(big_integer const &rhs) {
    if (rhs.is_zero()) {
        return *this = ZERO;
    }
    const big_integer copy = *this;
    data = std::vector<ui>(1, 0);
    big_integer temp;
    temp.len_prepare(rhs.length() + 1);
    for (size_t j = 0; j < rhs.length(); j++) {
        temp.mul_short(copy, rhs.get_digit(j));
        add_abs(temp, j);
    }
    sign = copy.sign ^ rhs.sign;
    shrink_zeroes();
    return *this;
}

ui big_integer::div_short(const unsigned int &a) {
    if (!a) {
        throw "Division by zero!";
    }
    ui64 cf = 0, res;
    for (size_t i = data.size(); i > 0; i--) {
        res = suic(data[i - 1]) + (cf << POWER);
        data[i - 1] = suic(res / a);
        cf = (res % a);
    }
    shrink_zeroes();
    return suic(cf);
}

big_integer big_integer::divide(const big_integer &rhs, const bool &mode) {
    if (rhs.is_zero()) {
        throw std::runtime_error("Division by zero");
    }

    bool ls = sign, rs = rhs.sign;
    sign = false;
    int8_t cr = comp_abs(*this, rhs);
    if (cr == 0) {
        if (mode) {
            return *this;
        } else {
            return ONE;
        }
    }
    if (cr == -1) {
        if (mode) {
            return *this;
        } else {
            return ZERO;
        }
    }

    if (rhs.length() == 1) {
        ui rem = div_short(rhs.get_digit(0));
        sign = ls ^ rs;
        if (mode) {
            big_integer rem_big(rem);
            rem_big.sign = ls;
            return rem_big;
        } else {
            return *this;
        }
    }

    big_integer b(rhs);
    b.sign = false;
    ui d;
    if (b.data.back() < HALFBASE && b.data.back() > HALFBASE / 2) {
        d = 1;
    } else if (b.data.back() > HALFBASE) {
        d = 0;
    } else {
        d = suic(ceil(log2(static_cast<double>(HALFBASE) / b.data.back())));
    }
    b <<= d;
    *this <<= d;
    size_t m = length() - b.length();
    size_t n = b.length();
    big_integer q;
    q.data.resize(m + 1, 0);
    big_integer r;
    r.data.resize(n, 0);

    for (size_t i = 0; i < n; i++) {
        r.data[i] = data[length() - n + i];
    }

    if (r > b) {
        q.data[m] = 1;
        r -= b;
    }

    big_integer c;
    size_t j;
    for (size_t i = m; i > 0; i--) {
        j = i - 1;
        r <<= POWER;
        r.data[0] = data[length() - n + j - m];
        q.data[j] = suic(((sui64c(r.get_digit(r.length() - 1)) << POWER) + r.get_digit(r.length() - 2)) /
                         sui64c(b.data.back()));

        c = b;
        c.mul_short(c, q.data[j]);
        r -= c;
        while (r.sign) {
            q.data[j]--;
            r += b;
        }
    }

    q.shrink_zeroes();
    r.shrink_zeroes();
    r.sign = ls;
    q.sign = ls ^ rs;
    if (mode) {
        return r >> d;
    } else {
        return q;
    }
}

big_integer &big_integer::operator/=(big_integer const &rhs) {
    return *this = divide(rhs, false);
}

big_integer &big_integer::operator%=(big_integer const &rhs) {
    return *this = divide(rhs, true);
}

void big_integer::bit_prepare(const big_integer &r) {
    len_prepare(r.length());
    if (sign) {
        complem_two();
    }
}

template<class FunctorT>
big_integer &big_integer::bit_gen(const big_integer &rhs, FunctorT functor) {

    big_integer r = rhs;
    r.bit_prepare(*this);
    bit_prepare(r);
    for (size_t i = 0; i < rhs.length(); i++) {
        data[i] = functor(get_digit(i), r.get_digit(i));
    }
    ~r;
    sign = functor(sign, rhs.sign);
    if (sign) {
        recomplem_two();
    }
    shrink_zeroes();
    return *this;
}

big_integer &big_integer::operator&=(big_integer const &rhs) {
    return bit_gen(rhs, std::bit_and<ui>());
}

big_integer &big_integer::operator|=(big_integer const &rhs) {
    return bit_gen(rhs, std::bit_or<ui>());
}

big_integer &big_integer::operator^=(big_integer const &rhs) {
    return bit_gen(rhs, std::bit_xor<ui>());
}

big_integer &big_integer::operator<<=(int rhs) {
    size_t fulls = rhs / POWER;
    size_t shift = rhs % POWER;
    ui64 cf = 0;
    size_t os = data.size();
    data.resize(os + fulls + 1, 0);
    for (size_t i = os; i > 0; i--) {
        cf = sui64c(data[i + fulls]) << POWER;
        cf |= sui64c(data[i - 1]) << shift;
        data[i + fulls] = suic(cf >> POWER);
        data[i - 1 + fulls] = suic(cf);
    }
    for (size_t i = fulls; i > 0; i--) {
        data[i - 1] = 0;
    }
    shrink_zeroes();
    return *this;
}

big_integer &big_integer::operator>>=(int rhs) {
    size_t fulls = rhs / POWER;
    size_t shift = rhs % POWER;
    ui64 cf = 0;

    for (size_t i = 0; i < length() - fulls; i++) {
        cf = sui64c(data[i]);
        cf |= sui64c(get_digit(i + 1)) << POWER;
        data[i] = suic(cf >> shift);
    }

    data.resize(length() - fulls);

    shrink_zeroes();
    if (sign) {
        *this -= 1;
    }
    return *this;
}

big_integer big_integer::operator+() const {
    return *this;
}

big_integer big_integer::operator-() const {
    big_integer r = *this;
    r.sign ^= 1;
    return r;
}

big_integer big_integer::operator~() const {
    big_integer r = *this;
    r.sign ^= 1;
    r -= 1;
    return r;
}

big_integer &big_integer::operator++() {
    *this += 1;
    return *this;
}

big_integer big_integer::operator++(int) {
    big_integer r = *this;
    ++*this;
    return r;
}

big_integer &big_integer::operator--() {
    *this -= 1;
    return *this;
}

big_integer big_integer::operator--(int) {
    big_integer r = *this;
    --*this;
    return r;
}

big_integer operator+(big_integer a, big_integer const &b) {
    return a += b;
}

big_integer operator-(big_integer a, big_integer const &b) {
    return a -= b;
}

big_integer operator*(big_integer a, big_integer const &b) {
    return a *= b;
}

big_integer operator/(big_integer a, big_integer const &b) {
    return a /= b;
}

big_integer operator%(big_integer a, big_integer const &b) {
    return a %= b;
}

big_integer operator&(big_integer a, big_integer const &b) {
    return a &= b;
}

big_integer operator|(big_integer a, big_integer const &b) {
    return a |= b;
}

big_integer operator^(big_integer a, big_integer const &b) {
    return a ^= b;
}

big_integer operator<<(big_integer a, int b) {
    return a <<= b;
}

big_integer operator>>(big_integer a, int b) {
    return a >>= b;
}

bool operator==(big_integer const &a, big_integer const &b) {
    return comp(a, b) == 0;
}

bool operator!=(big_integer const &a, big_integer const &b) {
    return comp(a, b) != 0;
}

bool operator<(big_integer const &a, big_integer const &b) {
    return comp(a, b) < 0;
}

bool operator>(big_integer const &a, big_integer const &b) {
    return comp(a, b) > 0;
}

bool operator<=(big_integer const &a, big_integer const &b) {
    return comp(a, b) <= 0;
}

bool operator>=(big_integer const &a, big_integer const &b) {
    return comp(a, b) >= 0;
}

std::string to_string(big_integer const &a) {
    if (a.data.empty() || a.is_zero()) {
        return "0";
    }

    ui cur_digit;
    big_integer v = a;
    v.shrink_zeroes();
    std::string res;
    bool check;
    while (!v.is_zero()) {
        cur_digit = v.div_short(BASE_10);
        check = v.is_zero();
        for (int i = 0; i < 9 && (!check || cur_digit); i++) {
            res.push_back('0' + char(cur_digit % 10));
            cur_digit /= 10;
        }
    }
    if (a.sign)
        res.push_back('-');
    std::reverse(res.begin(), res.end());
    return res;

}

std::ostream &operator<<(std::ostream &s, big_integer const &a) {
    return s << to_string(a);
}
