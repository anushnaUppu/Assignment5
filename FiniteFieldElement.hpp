

namespace Cryptography
{

        namespace   detail
        {

            int EGCD(int a, int b, int& u, int &v)
            {
                u = 1;
                v = 0;
                int g = a;
                int u1 = 0;
                int v1 = 1;
                int g1 = b;
                while (g1 != 0)
                {
                    int q = g/g1;
                    int t1 = u - q*u1;
                    int t2 = v - q*v1;
                    int t3 = g - q*g1;
                    u = u1; v = v1; g = g1;
                    u1 = t1; v1 = t2; g1 = t3;
                }

                return g;
            }

            int InvMod(int x, int n)
            {
                x = x % n;
                int u,v,g,z;
                g = EGCD(x, n, u,v);
                if (g != 1)
                {

                    z = 0;
                }
                else
                {
                    z = u % n;
                }
                return z;
            }
        }


        template<int P>
        class   FiniteFieldElement
        {
            int     i_;

            void    assign(int i)
            {
                i_ = i;
                if ( i<0 )
                {

                    i_ = (i%P) + 2*P;
                }

                i_ %= P;
            }

            public:

                FiniteFieldElement()
                 : i_(0)
                {}

                explicit FiniteFieldElement(int i)
                {
                    assign(i);
                }

                FiniteFieldElement(const FiniteFieldElement<P>& rhs)
                 : i_(rhs.i_)
                {
                }


                int i() const { return i_; }

                FiniteFieldElement  operator-() const
                {
                    return FiniteFieldElement(-i_);
                }

                FiniteFieldElement& operator=(int i)
                {
                    assign(i);
                    return *this;
                }

                FiniteFieldElement<P>& operator=(const FiniteFieldElement<P>& rhs)
                {
                    i_ = rhs.i_;
                    return *this;
                }

                FiniteFieldElement<P>& operator*=(const FiniteFieldElement<P>& rhs)
                {
                    i_ = (i_*rhs.i_) % P;
                    return *this;
                }

                friend bool    operator==(const FiniteFieldElement<P>& lhs, const FiniteFieldElement<P>& rhs)
                {
                    return (lhs.i_ == rhs.i_);
                }

                friend bool    operator==(const FiniteFieldElement<P>& lhs, int rhs)
                {
                    return (lhs.i_ == rhs);
                }

                friend bool    operator!=(const FiniteFieldElement<P>& lhs, int rhs)
                {
                    return (lhs.i_ != rhs);
                }

                friend FiniteFieldElement<P> operator/(const FiniteFieldElement<P>& lhs, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( lhs.i_ * detail::InvMod(rhs.i_,P));
                }

                friend FiniteFieldElement<P> operator+(const FiniteFieldElement<P>& lhs, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( lhs.i_ + rhs.i_);
                }
                // a - b
                friend FiniteFieldElement<P> operator-(const FiniteFieldElement<P>& lhs, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( lhs.i_ - rhs.i_);
                }
                // a + int
                friend FiniteFieldElement<P> operator+(const FiniteFieldElement<P>& lhs, int i)
                {
                    return FiniteFieldElement<P>( lhs.i_+i);
                }
                // int + a
                friend FiniteFieldElement<P> operator+(int i, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( rhs.i_+i);
                }
                // int * a
                friend FiniteFieldElement<P> operator*(int n, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( n*rhs.i_);
                }
                // a * b
                friend FiniteFieldElement<P> operator*(const FiniteFieldElement<P>& lhs, const FiniteFieldElement<P>& rhs)
                {
                    return FiniteFieldElement<P>( lhs.i_ * rhs.i_);
                }

                template<int T>
                friend  ostream&    operator<<(ostream& os, const FiniteFieldElement<T>& g)
                {
                    return os << g.i_;
                }
        };
}


