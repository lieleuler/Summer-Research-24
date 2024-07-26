
classdef MCF
    properties (Access = private)
        terms
        max_terms
    end
    methods
        % == Self Logic == %
        function result = normalize(this)
            normalized_terms = this.terms;
            while length(normalized_terms) >= this.max_terms
                N = length(normalized_terms) - 1;
                s = normalized_terms(N);
                t = zeros(1, N);
                for i = N:-1:1
                    curr_sum = sum_floats_to_MCF(normalized_terms(i), s);
                    s = curr_sum(1);
                    t(i) = curr_sum(2);
                end
            
                k = 1;
                normalized_terms = zeros(1, N);
                for i = 1:N
                    curr_sum = sum_floats_to_MCF(s, t(i));
                    if curr_sum(2) ~= 0
                        normalized_terms(k) = curr_sum(1);
                        s = curr_sum(2);
                        k = k + 1;
                    end
                end
            end
            result = MCF(normalized_terms);
        end

        % == Constructor == %
        function my_obj = MCF(float)
            if nargin == 1
                my_obj.terms = float;
                my_obj.max_terms = 12;
            end
        end

        % == Arithmetic == %
        function result = plus(obj1, obj2)
            if isa(obj1, "MCF")
                obj1_terms = obj1.terms;
                max_terms = obj1.max_terms;
            else
                obj1_terms = obj1;
            end
            if isa(obj2, "MCF")
                obj2_terms = obj2.terms;
                max_terms = obj2.max_terms;
            else
                obj2_terms = obj2;
            end

            result = sum_MCF(obj1_terms, obj2_terms, max_terms);
        end

        function result = minus(obj1, obj2)
            result = obj1 + -obj2;
        end

        function result = mtimes(obj1, obj2)
            if isa(obj1, "MCF")
                obj1_terms = obj1.terms;
                max_terms = obj1.max_terms;
            else
                obj1_terms = obj1;
            end
            if isa(obj2, "MCF")
                obj2_terms = obj2.terms;
                max_terms = obj2.max_terms;
            else
                obj2_terms = obj2;
            end

            result = multiply_MCF(obj1_terms, obj2_terms, max_terms);
        end

        function result = mrdivide(obj1, obj2)
            result = double(obj1)/double(obj2);
        end

        function result = power(obj1, obj2)
            result = obj1;
            for i = 1:(obj2 - 1)
                result = multiply_MCF(result, obj1, obj1.max_terms);
            end
        end 

        function result = uminus(obj)
            result = MCF(-obj.terms);
        end

        % == Advanced Mathematical == %
        function result = sqrt(obj1)
            result = sqrt(double(obj1));
        end

        function result = abs(obj1)
            result = abs(double(obj1));
        end

        % == Logic == %
        function result = eq(obj1, obj2)
            result = double(obj1) == double(obj2);
        end

        function result = ge(obj1, obj2)
            result = double(obj1) >= double(obj2);
        end

        function result = le(obj1, obj2)
            result = double(obj1) <= double(obj2);
        end

        function result = gt(obj1, obj2)
            result = double(obj1) > double(obj2);
        end

        function result = lt(obj1, obj2)
            result = double(obj1) <= double(obj2);
        end

        % == Conversion == %
        function result = double(obj)
            result = sum(obj.terms);
        end

        function result = real(obj)
            result = MCF(real(obj.terms));
        end

        function result = imag(obj)
            result = MCF(imag(obj.terms));
        end
    end
end    