#include <bits/stdc++.h>
using namespace std;

struct Term {
    long long a; // coefficient
    int b, c, d; // exponents for x, sin, cos
    Term(long long a_=0,int b_=0,int c_=0,int d_=0):a(a_),b(b_),c(c_),d(d_){}
    bool sameKey(const Term& o) const {return b==o.b && c==o.c && d==o.d;}
};

struct Poly {
    vector<Term> terms;
    Poly(){}
    Poly(const Term& t){ if(t.a!=0) terms.push_back(t); }
    void addTerm(const Term& t){ if(t.a!=0) terms.push_back(t); }
    // simplify combine like terms
    void simplify(){
        sort(terms.begin(), terms.end(), [](const Term& x, const Term& y){
            if(x.b!=y.b) return x.b>y.b;
            if(x.c!=y.c) return x.c>y.c;
            return x.d>y.d;
        });
        vector<Term> res;
        for(auto &t: terms){
            if(!res.empty() && res.back().sameKey(t)){
                res.back().a += t.a;
            }else res.push_back(t);
        }
        // remove zero coeff
        terms.clear();
        for(auto &t: res) if(t.a!=0) terms.push_back(t);
    }
    // addition
    Poly operator+(const Poly& o) const{
        Poly r; r.terms = terms; r.terms.insert(r.terms.end(), o.terms.begin(), o.terms.end()); r.simplify(); return r; }
    // subtraction
    Poly operator-(const Poly& o) const{
        Poly r; r.terms = terms; for(auto &t: o.terms){ r.terms.emplace_back(-t.a, t.b, t.c, t.d); } r.simplify(); return r; }
    // multiplication
    Poly operator*(const Poly& o) const{
        Poly r; r.terms.reserve(terms.size()*o.terms.size());
        for(const auto &t1: terms){
            for(const auto &t2: o.terms){
                long long na = t1.a * t2.a;
                int nb = t1.b + t2.b;
                int nc = t1.c + t2.c;
                int nd = t1.d + t2.d;
                r.terms.emplace_back(na, nb, nc, nd);
            }
        }
        r.simplify();
        return r;
    }
    // derivative of term
    static Poly derivativeTerm(const Term& t){
        Poly res;
        // derivative of a*x^b*sin^c*cos^d
        // a * ( b*x^{b-1}*sin^c*cos^d )
        if(t.b!=0){
            res.addTerm(Term(t.a * t.b, t.b-1, t.c, t.d));
        }
        // a * ( c*x^b*sin^{c-1}*cos^{d+1} )
        if(t.c!=0){
            res.addTerm(Term(t.a * t.c, t.b, t.c-1, t.d+1));
        }
        // a * ( -d*x^b*sin^{c+1}*cos^{d-1} )
        if(t.d!=0){
            res.addTerm(Term(-t.a * t.d, t.b, t.c+1, t.d-1));
        }
        res.simplify();
        return res;
    }
    Poly derivative() const{
        Poly res;
        for(const auto &t: terms){
            Poly d = derivativeTerm(t);
            res = res + d;
        }
        res.simplify();
        return res;
    }
    // output string for poly (no outer parentheses)
    string toString() const{
        if(terms.empty()) return "0";
        string out;
        for(size_t i=0;i<terms.size();++i){
            const Term& t = terms[i];
            long long coeff = t.a;
            // sign handling
            if(i==0){
                if(coeff<0) out += "-"; // leading minus
            }else{
                out += (coeff<0?"-":"+");
            }
            long long abscoeff = llabs(coeff);
            bool needCoeff = (abscoeff!=1) || (t.b==0 && t.c==0 && t.d==0);
            if(needCoeff) out += to_string(abscoeff);
            // variables
            if(t.b!=0){ out += "x"; if(t.b!=1) out += "^" + to_string(t.b); }
            if(t.c!=0){ out += "sin"; if(t.c!=1) out += "^" + to_string(t.c); out += "x"; }
            if(t.d!=0){ out += "cos"; if(t.d!=1) out += "^" + to_string(t.d); out += "x"; }
        }
        return out;
    }
};

struct Frac {
    Poly p, q; // numerator / denominator
    Frac(){}
    Frac(const Poly& p_, const Poly& q_):p(p_),q(q_){}
    // constructors for constant
    static Frac fromConst(long long v){
        Poly num(Term(v,0,0,0));
        Poly den(Term(1,0,0,0));
        return Frac(num, den);
    }
    // addition
    Frac operator+(const Frac& o) const{
        Poly np = p * o.q + o.p * q;
        Poly nq = q * o.q;
        return Frac(np, nq);
    }
    // subtraction
    Frac operator-(const Frac& o) const{
        Poly np = p * o.q - o.p * q;
        Poly nq = q * o.q;
        return Frac(np, nq);
    }
    // multiplication
    Frac operator*(const Frac& o) const{
        return Frac(p * o.p, q * o.q);
    }
    // division
    Frac operator/(const Frac& o) const{
        return Frac(p * o.q, q * o.p);
    }
    // derivative
    Frac derivative() const{
        // (p' * q - p * q') / (q * q)
        Poly np = p.derivative() * q - p * q.derivative();
        Poly nq = q * q;
        return Frac(np, nq);
    }
    // simplify numerator and denominator (combine like terms)
    void simplify(){ p.simplify(); q.simplify(); }
    string toString() const{
        // handle special cases
        // if numerator is zero
        if(p.terms.empty()) return "0";
        // if denominator is 1
        bool denomIsOne = (q.terms.size()==1 && q.terms[0].a==1 && q.terms[0].b==0 && q.terms[0].c==0 && q.terms[0].d==0);
        string numStr = p.toString();
        string denStr = q.toString();
        // need parentheses if more than one term
        bool numMulti = p.terms.size()>1;
        bool denMulti = q.terms.size()>1;
        if(denomIsOne) return numStr;
        string res;
        if(numMulti) res += "(" + numStr + ")"; else res += numStr;
        res += "/";
        if(denMulti) res += "(" + denStr + ")"; else res += denStr;
        return res;
    }
};

// Tokenizer
enum TokenType {TOK_END, TOK_NUM, TOK_X, TOK_SIN, TOK_COS, TOK_PLUS, TOK_MINUS, TOK_MUL, TOK_DIV, TOK_LP, TOK_RP, TOK_POW};
struct Token{TokenType type; long long value;};

string s; size_t pos;
Token cur;

void skipSpaces(){ while(pos<s.size() && isspace(s[pos])) ++pos; }

bool startsWith(const string& pre){ return s.compare(pos, pre.size(), pre)==0; }

Token nextToken(){
    skipSpaces();
    if(pos>=s.size()) return {TOK_END,0};
    char c=s[pos];
    if(isdigit(c)){
        long long val=0; while(pos<s.size() && isdigit(s[pos])){ val=val*10+(s[pos]-'0'); ++pos; }
        return {TOK_NUM,val};
    }
    if(c=='x') { ++pos; return {TOK_X,0}; }
    if(startsWith("sin")) { pos+=3; return {TOK_SIN,0}; }
    if(startsWith("cos")) { pos+=3; return {TOK_COS,0}; }
    ++pos;
    switch(c){
        case '+': return {TOK_PLUS,0};
        case '-': return {TOK_MINUS,0};
        case '*': return {TOK_MUL,0};
        case '/': return {TOK_DIV,0};
        case '(': return {TOK_LP,0};
        case ')': return {TOK_RP,0};
        case '^': return {TOK_POW,0};
    }
    // unknown, skip
    return nextToken();
}

void advance(){ cur = nextToken(); }

// Forward declarations
Frac parseExpr();
Frac parseTerm();
Frac parseFactor();

int parseNumber(){ if(cur.type==TOK_NUM){ int v=cur.value; advance(); return v; } return 1; }

// parse primary which may be number, variable part, or parenthesized expression
Frac parsePrimary(){
    // unary minus
    if(cur.type==TOK_MINUS){
        advance();
        Frac val = parsePrimary();
        return Frac::fromConst(-1) * val;
    }
    if(cur.type==TOK_NUM){
        int val = cur.value; advance();
        return Frac::fromConst(val);
    }
    if(cur.type==TOK_X || cur.type==TOK_SIN || cur.type==TOK_COS){
        // implicit multiplication of possible sequence
        Poly resultPoly; bool first=true;
        while(cur.type==TOK_X || cur.type==TOK_SIN || cur.type==TOK_COS){
            Term t(1,0,0,0);
            if(cur.type==TOK_X){
                advance();
                int exp=1;
                if(cur.type==TOK_POW){ advance(); exp = parseNumber(); }
                t.b = exp;
            }else if(cur.type==TOK_SIN){
                advance();
                int exp=1;
                if(cur.type==TOK_POW){ advance(); exp = parseNumber(); }
                // consume the following 'x' as argument of sin
                if(cur.type==TOK_X) advance();
                t.c = exp;
            }else if(cur.type==TOK_COS){
                advance();
                int exp=1;
                if(cur.type==TOK_POW){ advance(); exp = parseNumber(); }
                // consume the following 'x' as argument of cos
                if(cur.type==TOK_X) advance();
                t.d = exp;
            }
            Poly p(t);
            if(first){ resultPoly = p; first=false; }
            else resultPoly = resultPoly * p;
        }
        return Frac(resultPoly, Poly(Term(1,0,0,0)));
    }
    if(cur.type==TOK_LP){ advance(); Frac inside = parseExpr(); if(cur.type==TOK_RP) advance(); return inside; }
    // fallback zero
    return Frac::fromConst(0);
}

// parse factor handles division and multiplication with explicit operators and also implicit multiplication of consecutive primaries (handled in parsePrimary)
Frac parseFactor(){
    Frac left = parsePrimary();
    while(true){
        if(cur.type==TOK_MUL){ advance(); Frac right = parsePrimary(); left = left * right; }
        else if(cur.type==TOK_DIV){ advance(); Frac right = parsePrimary(); left = left / right; }
        else if(cur.type==TOK_NUM || cur.type==TOK_X || cur.type==TOK_SIN || cur.type==TOK_COS || cur.type==TOK_LP){
            // implicit multiplication (e.g., 3x, 2sinx, (a)(b))
            Frac right = parsePrimary();
            left = left * right;
        }
        else break;
    }
    return left;
}

Frac parseTerm(){ // handles + -
    Frac left = parseFactor();
    while(true){
        if(cur.type==TOK_PLUS){ advance(); Frac right = parseFactor(); left = left + right; }
        else if(cur.type==TOK_MINUS){ advance(); Frac right = parseFactor(); left = left - right; }
        else break;
    }
    return left;
}

Frac parseExpr(){ return parseTerm(); }

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string input;
    if(!std::getline(cin, input)) return 0;
    // remove possible carriage return
    if(!input.empty() && input.back()=='\r') input.pop_back();
    // remove spaces from input
    input.erase(remove(input.begin(), input.end(), ' '), input.end());
    s = input; pos=0; advance();
    if(cur.type==TOK_PLUS) advance(); // skip leading plus
    Frac expr = parseExpr(); expr.simplify();
    Frac deriv = expr.derivative(); deriv.simplify();
    cout << expr.toString() << "\n" << deriv.toString() << "\n";
    return 0;
}
