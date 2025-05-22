module iso_varying_string
    !! This module defines facilities in Fortran for the manipulation of
    !! character strings of dynamically variable length. It is in conformance
    !! with the [ISO/IEC 1539-2: 2000](http://www.astro.wisc.edu/~townsend/resource/download/code/Fortran-ISO_VARYING_STRING.pdf)
    !! extension to the Fortran Standard.
    !!
    !! Neither the internal representation of the derived type, nor the
    !! algorithms used to implement the procedures or operators are directed
    !! by the standard, and as such should not be relied upon by any user of
    !! this module.
    !!
    !! It should be noted that this module defines facilities for dynamically
    !! varying length strings of characters of default kind only. Similar
    !! facilities could be defined for non-default kind characters by a
    !! separate, if similar, module for each such character kind.
    !!
    !! This module has been designed, as far as is reasonable, to provide for
    !! varying length character strings the facilities that are available for
    !! intrinsic fixed length character strings. All the intrinsic operations
    !! and functions that apply to fixed length character strings have extended
    !! meanings defined by this module for varying length character strings.
    !! Also, a small number facilities are defined that are appropriate because
    !! of the essential differences between the intrinsic type and the varying
    !! length derived data type.

    USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR

    implicit none
    private
    public :: &
            varying_string, &
            assignment(=), &
            operator(//), &
            operator(==), &
            operator(/=), &
            operator(<), &
            operator(<=), &
            operator(>), &
            operator(>=), &
            adjustl, &
            adjustr, &
            char, &
            iachar, &
            ichar, &
            index, &
            len, &
            len_trim, &
            lge, &
            lgt, &
            lle, &
            llt, &
            repeat, &
            scan, &
            trim, &
            verify, &
            var_str, &
            get, &
            put, &
            put_line, &
            extract, &
            insert, &
            remove, &
            replace, &
            split

    type :: varying_string ! Sec. 3.2
        !! Entities of this type shall represent values that are strings of
        !! characters of default kind. These character strings may be of any
        !! non-negative length and this length may vary dynamically during the
        !! execution of a program. There shall be no arbitrary upper length
        !! limit other than that imposed by the size of the processor and the
        !! complexity of the programs it is able to process. The characters
        !! representing the value of the string have positions 1,2,...,N, where
        !! N is the length of the string. The internal structure of the type
        !! shall be `PRIVATE` to the module.
        private
        character(len=1), allocatable :: characters(:)
    contains
        private
        procedure :: write_formatted
        procedure :: write_unformatted
        procedure :: read_formatted
        procedure :: read_unformatted
        generic, public :: write(formatted) => write_formatted
            !! Enables use in formatted (list-directed, namelist and explicit
            !! format) write and print statements
        generic, public :: write(unformatted) => write_unformatted
            !! Enables use in unformatted write statements
        generic, public :: read(formatted) => read_formatted
            !! Enables use in formatted (list-directed, namelist and explicit
            !! format) read statements
        generic, public :: read(unformatted) => read_unformatted
            !! Enable use in unformatted read statements
    end type

    interface assignment(=) ! Sec. 3.3.1
        !! An elemental assignment of the form
        !!
        !! ```fortran
        !! var = expr
        !! ```
        !!
        !! shall be defined with following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! **Action**. The characters that are the value of the expression
        !! `expr` become the value of the variable `var`. There are two cases:
        !!
        !! * *Case(i):* Where the variable is of type `VARYING_STRING`, the
        !!   the length of the variable becomes that of the expression.
        !! * *Case(i):* Where the variable is of type `CHARACTER`, the rules
        !!   of intrinsic assignment to a Fortran character variable apply.
        !!   Namely, if the expression string is longer than the declared
        !!   length of the character variable, only the left-most characters are
        !!   assigned. If the character variable is longer than that of the
        !!   string expression, it is padded on the right with blanks.
        module procedure assign_character_to_string
        module procedure assign_string_to_character
    end interface

    interface operator(//) ! Sec. 3.3.2
        !! The elemental concatenation operation
        !!
        !! ```fortran
        !! string_a // string_b
        !! ```
        !!
        !! shall be defined with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is a new string whose characters
        !! are the same as those produced by concatenating the operand strings
        !! in the order given.
        module procedure concat_strings
        module procedure concat_string_and_character
        module procedure concat_character_and_string
    end interface

    interface operator(==) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a == string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.EQ.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` contains
        !! the same characters in the same sequence as `string_b` and false
        !! otherwise. If `string_a` and `string_b` are of different lengths, the
        !! comparison is done as if the shorter string were padded on the right
        !! with blanks.
        module procedure string_eq_string
        module procedure character_eq_string
        module procedure string_eq_character
    end interface

    interface operator(/=) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a /= string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.NE.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is false if `string_a` contains
        !! the same characters in the same sequence as `string_b` and true
        !! otherwise. If `string_a` and `string_b` are of different lengths, the
        !! comparison is done as if the shorter string were padded on the right
        !! with blanks.
        module procedure string_ne_string
        module procedure character_ne_string
        module procedure string_ne_character
    end interface

    interface operator(<) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a < string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.LT.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is
        !! lexically less than `string_b` and false otherwise. The collating
        !! sequence used for the comparisons is that defined by the processor
        !! for characters of default kind. If `string_a` and `string_b` are of
        !! different lengths, the comparison is done as if the shorter string
        !! were padded on the right with blanks.
        module procedure string_lt_string
        module procedure character_lt_string
        module procedure string_lt_character
    end interface

    interface operator(<=) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a <= string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.LT.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is
        !! lexically less than or equal to `string_b` and false otherwise. The
        !! collating sequence used for the comparisons is that defined by the
        !! processor for characters of default kind. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks.
        module procedure string_le_string
        module procedure character_le_string
        module procedure string_le_character
    end interface

    interface operator(>) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a > string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.LT.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is
        !! lexically greater than `string_b` and false otherwise. The collating
        !! sequence used for the comparisons is that defined by the processor
        !! for characters of default kind. If `string_a` and `string_b` are of
        !! different lengths, the comparison is done as if the shorter string
        !! were padded on the right with blanks.
        module procedure string_gt_string
        module procedure character_gt_string
        module procedure string_gt_character
    end interface

    interface operator(>=) ! Sec. 3.3.3
        !! Elemental comparison of the form
        !!
        !! ```fortran
        !! string_a >= string_b
        !! ```
        !!
        !! shall be defined for operands with the following type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! The values of the operands are unchanged by the operation. Note that
        !! the equivalent operator form, `.LT.`, also has its meaning extended
        !! in this manner.
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is
        !! lexically greater than or equal to `string_b` and false otherwise.
        !! The collating sequence used for the comparisons is that defined by
        !! the processor for characters of default kind. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks.
        module procedure string_ge_string
        module procedure character_ge_string
        module procedure string_ge_character
    end interface

    interface adjustl ! Sec. 3.4.1
        !! **Description**. Adjusts to the left, removing any leading blanks and
        !! inserting trailing blanks.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `string` shall be of type `VARYING_STRING`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is the same as `string`, except
        !! that any leading blanks have been deleted and the same number of
        !! trailing blanks inserted.
        module procedure string_adjustl
    end interface

    interface adjustr ! Sec. 3.4.1
        !! **Description**. Adjusts to the right, removing any trailing blanks
        !! and inserting leading blanks.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `string` shall be of type `VARYING_STRING`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is the same as `string`, except
        !! that any trailing blanks have been deleted and the same number of
        !! leading blanks inserted.
        module procedure string_adjustr
    end interface

    interface char ! Sec. 3.4.3
        !! **Description.** Converts a varying string value to default `CHARACTER`.
        !!
        !! **Class.** Pure transformational function.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be scalar and of type `VARYING_STRING`.
        !! * `length` (optional) shall be scalar and of type default `INTEGER`.
        !!
        !! **Result Characteristics.** Scalar of type default `CHARACTER`. If
        !! `length` is absent, the result has the same length as `string`. If
        !! `length` is present, the result has the length specified by the `length`.
        !!
        !! **Result Value.**
        !!
        !! * *Case(i):* If `length` is absent, the result is a copy of the
        !!   characters in the argument `string`.
        !! * *Case(ii):* If `length` is present, the result is a copy of the
        !!   characters in the argument `string`, that may have been truncated
        !!   or padded. If `string` is longer than `length`, the result is
        !!   truncated on the right. If `string` is shorter than `length`, the
        !!   result is padded on the right with blanks. If `length` is less than
        !!   one, the result is of zero length.
        module procedure string_to_char
        module procedure string_to_char_with_length
    end interface

    interface iachar ! Sec. 3.4.4
        !! **Description.** Returns the position of a character in the collating
        !! sequence defined by the International Standard ISO 646: 1991.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `c` shall be of type `VARYING_STRING` and of length
        !! exactly one.
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.** The result value is the position of the character
        !! `c` in the collating sequence defined by the International Standard
        !! ISO 646: 1991 for default characters. If the character `c` is not
        !! defined in the standard set, the result is processor dependent but is
        !! always equal to `IACHAR(CHAR(c))`.
        module procedure string_iachar
    end interface

    interface ichar ! Sec. 3.4.5
        !! **Description.** Returns the position of a character in the processor
        !! defined collating sequence.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `c` shall be of type `VARYING_STRING` and of length
        !! exactly one.
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.** The result value is the position of the character
        !! `c` in the processor defined collating sequence for default
        !! characters. That is the result value is `ICHAR(CHAR(c))`.
        module procedure string_ichar
    end interface

    interface index ! Sec. 3.4.6
        !! **Description.** Returns an integer that is the starting position of
        !! a substring within a string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! `string` and `substring` shall be of one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! `back` (optional) shall be of type default `LOGICAL`.
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.**
        !!
        !! * *Case(i):* If `back` is absent or is present with the value false,
        !!   the result is the minimum positive value of `I` such that
        !!   `EXTRACT(string,I,I+LEN(substring)-1)==substring`, (where `EXTRACT`
        !!   is defined in this module) or zero if there is no such value.
        !! * *Case(ii):* If `back` is present with the value true, the result is
        !!   the maximum value of `I` less than or equal to
        !!   `LEN(string)-LEN(substring)+1` such that
        !!   `EXTRACT(string,I,I+LEN(substring))==substring`, or zero if there
        !!   is no such value.
        module procedure string_index_string
        module procedure string_index_character
        module procedure character_index_string
    end interface

    interface len ! Sec. 3.4.7
        !! **Description.** Returns the length of a character string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `string` shall be of type `VARYING_STRING`.
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.** The result value is the number of characters in
        !! `string`.
        !!
        !! **Note.** This function is not elemental for `string` of type
        !! `CHARACTER`.
        module procedure len_string
    end interface

    interface len_trim ! Sec. 3.4.8
        !! **Description.** Returns the length of a string not counting any
        !! trailing blanks.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `string` shall be of type `VARYING_STRING`.
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.** The result value is the position of the last
        !! non-blank character in `string`. If the argument `string` contains
        !! only blank characters or is of zero length, the result is zero.
        module procedure len_trim_string
    end interface

    interface lge ! Sec. 3.4.9
        !! **Description.** Compares the lexical ordering of two strings based
        !! on the ISO 646: 1991 collating sequence.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string_a` and `string_b` shall be of one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is lexically
        !! greater than or equal to `string_b`, and is false otherwise. The
        !! collating sequence used to establish the ordering of characters is
        !! that of the International Standard ISO 646: 1991. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks. If either
        !! argument contains a character `c` not defined by the standard, the
        !! result value is processor dependent and based on the collating value
        !! for `IACHAR(c)`. Zero length strings are considered to be lexically
        !! equal.
        module procedure string_lge_string
        module procedure character_lge_string
        module procedure string_lge_character
    end interface

    interface lgt ! Sec. 3.4.10
        !! **Description.** Compares the lexical ordering of two strings based
        !! on the ISO 646: 1991 collating sequence.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string_a` and `string_b` shall be of one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is lexically
        !! greater than `string_b`, and is false otherwise. The
        !! collating sequence used to establish the ordering of characters is
        !! that of the International Standard ISO 646: 1991. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks. If either
        !! argument contains a character `c` not defined by the standard, the
        !! result value is processor dependent and based on the collating value
        !! for `IACHAR(c)`. Zero length strings are considered to be lexically
        !! equal.
        module procedure string_lgt_string
        module procedure character_lgt_string
        module procedure string_lgt_character
    end interface

    interface lle ! Sec. 3.4.11
        !! **Description.** Compares the lexical ordering of two strings based
        !! on the ISO 646: 1991 collating sequence.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string_a` and `string_b` shall be of one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is lexically
        !! less than or equal to `string_b`, and is false otherwise. The
        !! collating sequence used to establish the ordering of characters is
        !! that of the International Standard ISO 646: 1991. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks. If either
        !! argument contains a character `c` not defined by the standard, the
        !! result value is processor dependent and based on the collating value
        !! for `IACHAR(c)`. Zero length strings are considered to be lexically
        !! equal.
        module procedure string_lle_string
        module procedure character_lle_string
        module procedure string_lle_character
    end interface

    interface llt ! Sec. 3.4.12
        !! **Description.** Compares the lexical ordering of two strings based
        !! on the ISO 646: 1991 collating sequence.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string_a` and `string_b` shall be of one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! **Result Characteristics.** Of type default `LOGICAL`.
        !!
        !! **Result Value.** The result value is true if `string_a` is lexically
        !! less than `string_b`, and is false otherwise. The
        !! collating sequence used to establish the ordering of characters is
        !! that of the International Standard ISO 646: 1991. If `string_a` and
        !! `string_b` are of different lengths, the comparison is done as if the
        !! shorter string were padded on the right with blanks. If either
        !! argument contains a character `c` not defined by the standard, the
        !! result value is processor dependent and based on the collating value
        !! for `IACHAR(c)`. Zero length strings are considered to be lexically
        !! equal.
        module procedure string_llt_string
        module procedure character_llt_string
        module procedure string_llt_character
    end interface llt

    interface repeat ! Sec. 3.4.13
        !! **Description.** Concatenates several copies of a string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! * `string` - shall be of type `VARYING_STRING`.
        !! * `ncopies` - shall be of type default `INTEGER`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`
        !!
        !! **Result Value.** The result value is the string produced by repeated
        !! concatenation of the argument `string`, producing a string containing
        !! `ncopies` copies of `string`. If the value of `ncopies` is not
        !! positive, the result is of zero length.
        !!
        !! **Note.** This function is not elemental for `string` of type
        !! `CHARACTER`.
        module procedure string_repeat
    end interface

    interface scan ! Sec. 3.4.14
        !! **Description.** Scans a string for any one of the characters in a
        !! set of characters.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string` and `set` shall be one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! `back` (optional) shall be of type default `LOGICAL`
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.**
        !!
        !! * *Case(i):* If `back` is absent or is present with the value false
        !!   and if `string` contains at least one character that is in `set`,
        !!   the value of the result is the position of the left-most character
        !!   of `string` that is in `set`.
        !! * *Case(ii):* If `back` is present with the value true and if
        !!   `string` contains at least one character that is in `set`, the
        !!   value of the result is the position of the right-most character of
        !!   `string` that is in `set`.
        !! * *Case(iii):* The value of the result is zero if no charaqcter of
        !!   `string` is in `set` or if the length of either `string` or `set`
        !!   is zero.
        module procedure string_scan_string
        module procedure string_scan_character
        module procedure character_scan_string
    end interface

    interface trim ! Sec. 3.4.15
        !! **Description.** Removes trailing blanks from a string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `string` shall be of type `VARYING_STRING`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is the same as `string` except
        !! that any trailing blanks have been deleted. If the argument `string`
        !! contains only blank characters or is of zero length, the result is a
        !! zero-length string.
        !!
        !! **Note.** This function is not elemental for `string` of type
        !! `CHARACTER`.
        module procedure trim_string
    end interface

    interface verify ! Sec. 3.4.16
        !! **Description.** Verifies that a string contains only characters from
        !! a given set by scanning for any character not in the set.
        !!
        !! **Class.** Elemental function
        !!
        !! **Arguments.**
        !!
        !! `string` and `set` shall be one of the type combinations:
        !!
        !! * `VARYING_STRING` and `VARYING_STRING`
        !! * `VARYING_STRING` and `CHARACTER`
        !! * `CHARACTER` and `VARYING_STRING`
        !!
        !! `back` (optional) shall be of type default `LOGICAL`
        !!
        !! **Result Characteristics.** Of type default `INTEGER`.
        !!
        !! **Result Value.**
        !!
        !! * *Case(i):* If `back` is absent or is present with the value false
        !!   and if `string` contains at least one character that is not in
        !!   `set`, the value of the result is the position of the left-most
        !!   character of `string` that is not in `set`.
        !! * *Case(ii):* If `back` is present with the value true and if
        !!   `string` contains at least one character that is not in `set`, the
        !!   value of the result is the position of the right-most character of
        !!   `string` that is not in `set`.
        !! * *Case(iii):* The value of the result is zero if each character of
        !!   `string` is in `set` or if the length of `string` is zero.
        module procedure string_verify_string
        module procedure string_verify_character
        module procedure character_verify_string
    end interface

    interface var_str
        !! **Description.** Converts an intrinsic fixed-length character value
        !! into the equivalent varying-length string value.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Argument.** `char` shall be of type default `CHARACTER` and may be
        !! of any length.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is the same string of characters
        !! as the argument.
        module procedure var_str_char
    end interface

    interface get ! Sec. 3.6.1
        !! **Description.** Reads characters from an external file into a
        !! string.
        !!
        !! **Class.** Subroutine.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be scalar and of type `VARYING_STRING`. It is an
        !!   `INTENT(OUT)` argument.
        !! * `maxlen` (optional) shall be scalar and of type default `INTEGER`.
        !!   It is an `INTENT(IN)` argument.
        !! * `unit` shall be scalar and of type default `INTEGER`. It is an
        !!   `INTENT(IN)` argument that specifies the input unit to be used. The
        !!   unit shall be connected to a formatted file for sequential read
        !!   access. If the argument `unit` is omitted, the default input unit
        !!   is used.
        !! * `set` shall be scalar and either of type `VARYING_STRING` or of
        !!   type `CHARACTER`. It is an `INTENT(IN)` argument.
        !! * `separator` (optional) shall be scalar and of type
        !!   `VARYING_STRING`. It is an `INTENT(OUT)` argument.
        !! * `iostat` (optional) shall be scalar and of type default `INTEGER`.
        !!   It is an `INTENT(OUT)` argument.
        !!
        !! **Action.** The `GET` procedure causes characters from the connected
        !! file, starting with the next character in the current record if there
        !! is a current record or the first character of the next record if not,
        !! to be read and stored in the variable `string`. The end of record
        !! always terminates the input but may be terminated before this. If
        !! `maxlen` is present, its value indicates the maximum number of
        !! characters that will be read. If `maxlen` is less than or equal to
        !! zero, no characters will be read and `string` will be set to zero
        !! length. If `maxlen` is absent, a maximum of `HUGE(1)` is used. If the
        !! argument `set` is provided, this specifies a set of characters the
        !! occurrence of any of which will terminate the input. This terminal
        !! character, although read from the input file, will not be included in
        !! the result string. The file position after the data transfer is
        !! complete, is after the last character that was read. If the argument
        !! `separator` is present, the actual character found which terminates
        !! the transfer is returned in `separator`. If the transfer is
        !! terminated other than by the occurrence of a character in `set`, a
        !! zero length string is returned in `separator`. If the transfer is
        !! terminated by the end of record being reached, the file is positioned
        !! after the record just read. If present, the argument `iostat` is used
        !! to return the status resulting from the data transfer. A zero value
        !! is returned if a valid read operation occurs and the end-of-record is
        !! not reached, a positive value if an error occurs, and a negative
        !! value if an end-of-file or end-of-record condition occurs. Note, the
        !! negative value returned for an end-of-file condition shall be
        !! different from that returned for an end-of-record condition. If
        !! `iostat` is absent and an error or end-of-file condition occurs, the
        !! program execution is terminated.
        module procedure get_default_unit_to_end_of_record
        module procedure get_with_unit_to_end_of_record
        module procedure get_default_unit_to_terminator_string
        module procedure get_with_unit_to_terminator_string
        module procedure get_default_unit_to_terminator_characters
        module procedure get_with_unit_to_terminator_characters
    end interface

    interface put ! Sec. 3.6.2
        !! **Description.** Writes to an external file.
        !!
        !! **Class.** Subroutine.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be scalar and of type `VARYING_STRING` or type
        !!   `CHARACTER`. It is an `INTENT(IN)` argument.
        !! * `unit` shall be scalar and of type default `INTEGER`. It is an
        !!   `INTENT(IN)` argument that specifies the output unit to be used.
        !!   The unit shall be connected to a formatted file for sequential
        !!   write access. If the argument `unit` is omitted, the default output
        !!   unit is used.
        !! * `iostat` (optional) shall be scalar and of type default `INTEGER`.
        !!   It is an `INTENT(OUT)` argument.
        !!
        !! **Action.** The `PUT` procedure causes the characters of `string` to
        !! be appended to the current record, if there is a current record, or
        !! to the start of the next record if there is no current record. The
        !! last character transferred becomes the last character of the current
        !! record. If present, the argument `iostat` is used to return the
        !! status resulting from the data transfer. A zero value is returned if
        !! a valid write operation occurs, and a positive value if an error
        !! occurs. If `iostat` is absent and anything other than a valid write
        !! operation occurs, the program execution is terminated.
        module procedure put_String_Default_Unit
        module procedure put_string_with_unit
        module procedure put_characters_default_unit
        module procedure put_characters_with_unit
    end interface

    interface put_line ! Sec. 3.6.3
        !! **Description.** Writes to an external file and ends the record.
        !!
        !! **Class.** Subroutine.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be scalar and of type `VARYING_STRING` or type
        !!   `CHARACTER`. It is an `INTENT(IN)` argument.
        !! * `unit` shall be scalar and of type default `INTEGER`. It is an
        !!   `INTENT(IN)` argument that specifies the output unit to be used.
        !!   The unit shall be connected to a formatted file for sequential
        !!   write access. If the argument `unit` is omitted, the default output
        !!   unit is used.
        !! * `iostat` (optional) shall be scalar and of type default `INTEGER`.
        !!   It is an `INTENT(OUT)` argument.
        !!
        !! **Action.** The `PUT_LINE` procedure causes the characters of
        !! `string` to be appended to the current record, if there is a current
        !! record, or to the start of the next record if there is no current
        !! record. Following completion of the data transfer, the file is
        !! position after the record just written. If present, the argument
        !! `iostat` is used to return the status resulting from the data
        !! transfer. A zero value is returned if a valid write operation occurs,
        !! and a positive value if an error occurs. If `iostat` is absent and
        !! anything other than a valid write operation occurs, the program
        !! execution is terminated.
        module procedure put_line_string_default_unit
        module procedure put_line_string_with_unit
        module procedure put_line_characters_default_unit
        module procedure put_line_characters_with_unit
    end interface

    interface extract ! Sec. 3.7.1
        !! **Description.** Extracts a specified substring from a string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`.
        !! * `start` (optional) shall be of type default `INTEGER`.
        !! * `finish` (optional) shall be of type default `INTEGER`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is a copy of the characters of the
        !! argument `string` between positions `start` and `finish`, inclusive.
        !! If `start` is absent or less than one, the value one is used for
        !! `start`. If `finish` is absent or greater than `LEN(string)`, the
        !! value `LEN(string)` is used for `finish`. If `finish` is less than
        !! `start`, the result is a zero-length string.
        module procedure extract_character
        module procedure extract_string
    end interface

    interface insert ! Sec. 3.7.2
        !! **Description.** Inserts a substring into a string at a specified
        !! position.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`.
        !! * `start` shall be of type default `INTEGER`.
        !! * `substring` shall be either type `VARYING_STRING` or type default
        !!   `CHARACTER`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is a copy of the characters of the
        !! argument `string` with the characters of `substring` inserted into
        !! the copy of `string` before the character at the position `start`.
        !! If `start` is greater than `LEN(string)`, the value `LEN(string)+1`
        !! is used for `start` and substring is appended to the copy of
        !! `string`. If `start` is less than one, the value one is used for
        !! `start` and `substring` is inserted before the first character of
        !! the copy of `string`.
        module procedure insert_character_into_character
        module procedure insert_character_into_string
        module procedure insert_string_into_character
        module procedure insert_string_into_string
    end interface

    interface remove ! Sec. 3.7.3
        !! **Description.** Removes a specified substring from a string.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`.
        !! * `start` (optional) shall be of type default `INTEGER`.
        !! * `finish` (optional) shall be of type default `INTEGER`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is a copy of the characters of
        !! `string` with the characters between positions `start` and `finish`,
        !! inclusive, removed. If `start` is absent or less than one, the value
        !! one is used for `start`. If `finish` is absent or greater than
        !! `LEN(string)`, the value `LEN(string)` is used for `finish`. If
        !! `finish` is less than `start`, the characters of `string` are
        !! delivered unchanged as the result.
        module procedure remove_character
        module procedure remove_string
    end interface

    interface replace ! Sec. 3.7.4
        !! **Description.** Replaces a subset of the characters in a string by
        !! a given substring. The subset may be specified either by position or
        !! by content.
        !!
        !! **Class.** Elemental function.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`.
        !! * `start` shall be of type default `INTEGER`.
        !! * `finish` shall be of type default `INTEGER`.
        !! * `target` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`. It shall not be of zero length.
        !! * `every` (optional) shall be of type default `LOGICAL`.
        !! * `back` (optional) shall be of type default `LOGICAL`.
        !!
        !! **Result Characteristics.** Of type `VARYING_STRING`.
        !!
        !! **Result Value.** The result value is a copy of the characters in
        !! `string` modified as per one of the cases below.
        !!
        !! * *Case(i):* For a reference of the form
        !!   `REPLACE(string,start,substring)` the characters of the agrument
        !!   `substring` are inserted into the copy of `string` beginning with
        !!   the character at the character position `start`. The characters in
        !!   positions from `start` to `MIN(start+LEN(substring)-1,LEN(string))`
        !!   are deleted. If `start` is greater than `LEN(string)`, the value
        !!   `LEN(string)+1` is used for `start` and `substring` is appended to
        !!   the copy of `string`. If `start` is less than one, the value one
        !!   is used for `start`.
        !! * *Case(ii):* For a reference of the form
        !!   `REPLACE(string,start,finish,substring)` the characters in the copy
        !!   of `string` between positions `start` and `finish`, including those
        !!   at `start` and `finish`, are deleted and replaced by the
        !!   characters of `substring`. If `start` is less than one, the value
        !!   one is used for `start`. If `finish` is greater than `LEN(string)`,
        !!   the value `LEN(string)` is used for `finish`. If `finish` is less
        !!   than `start`, the characters of `substring` are inserted before
        !!   the character at `start` and no characters are deleted.
        !! * *Case(iii):* For a reference of the form
        !!   `REPLACE(string,target,substring,every,back)` the copy of `string`
        !!   is searched for occurrences of `target`. The search is done in the
        !!   `backward` direction if the argument `back` is present with the
        !!   value true, and in the forward direction otherwise. If `target` is
        !!   found, it is replaced by `substring`. If `every` is present with
        !!   the value true, the search and replace is continued from the
        !!   character following `target` in the search direction specified
        !!   until all occurrences of `target` in the copy of string are
        !!   replaced; otherwise only the first occurrence of `target` is
        !!   replaced.
        module procedure replace_character_with_character_start
        module procedure replace_string_with_character_start
        module procedure replace_character_with_string_start
        module procedure replace_string_with_string_start
        module procedure replace_character_with_character_range
        module procedure replace_string_with_character_range
        module procedure replace_character_with_string_range
        module procedure replace_string_with_string_range
        module procedure replace_target_character_with_character_in_character
        module procedure replace_target_character_with_character_in_string
        module procedure replace_target_character_with_string_in_character
        module procedure replace_target_character_with_string_in_string
        module procedure replace_target_string_with_character_in_character
        module procedure replace_target_string_with_character_in_string
        module procedure replace_target_string_with_string_in_character
        module procedure replace_target_string_with_string_in_string
    end interface

    interface split ! Sec. 3.7.5
        !! **Description.** Splits a string into two substrings with the
        !! substrings separated by the occurrence of a character from a
        !! specified separator set.
        !!
        !! **Class.** Elemental subroutine.
        !!
        !! **Arguments.**
        !!
        !! * `string` shall be type `VARYING_STRING`. It is an `INTENT(INOUT)`
        !!   argument.
        !! * `word` shall be of type `VARYING_STRING`. It is an `INTENT(OUT)`
        !!   argument.
        !! * `set` shall be either of type `VARYING_STRING` or type default
        !!   `CHARACTER`. It is an `INTENT(IN)` argument.
        !! * `separator` (optional) shall be of type `VARYING_STRING`. It is an
        !!   `INTENT(OUT)` argument.
        !! * `back` (optional) shall of type default `LOGICAL`. It is an
        !!   `INTENT(IN)` argument.
        !!
        !! **Action.** The effect of the procedure is to divide the `string` at
        !! the first occurrence of a character that is in `set`. The `string` is
        !! searched in the forward direction unless `back` is present with the
        !! value true, in which case the search is in the backward direction.
        !! The characters passed over in the search are returned in the argument
        !! `word` and the remainder of the string, not including the separator
        !! character, is returned in the argument `string`. If the argument
        !! `separator` is present, the actual character found which separates
        !! the `word` from the remainder of the `string` is returned in
        !! `separator`. If no character from `set` is found or `set` is of zero
        !! length, the whole string is returned in `word`, `string` is returned
        !! as zero length, and `separator` (if present) is returned as zero
        !! length. The effect of the procedure is such that, on return, either
        !! `word//separator//string` is the same as the initial string for a
        !! forward search, or `string//separator//word` is the same as the
        !! initial string for a backward search.
        module procedure split_character
        module procedure split_string
    end interface
contains
    elemental subroutine assign_character_to_string(lhs, rhs)
        ! Sec. 3.3.1
        type(varying_string), intent(out) :: lhs
        character(len=*), intent(in) :: rhs

        integer :: i
        integer :: length

        length = len(rhs)
        allocate(lhs%characters(length))
        do concurrent (i = 1 : length)
            lhs%characters(i) = rhs(i:i)
        end do
    end subroutine

    elemental subroutine assign_string_to_character(lhs, rhs)
        ! Sec. 3.3.1
        character(len=*), intent(out) :: lhs
        type(varying_string), intent(in) :: rhs

        integer :: i
        integer :: length_input
        integer :: length_output

        length_output = len(lhs)
        if (allocated(rhs%characters)) then
            length_input = size(rhs%characters)
            do concurrent (i = 1 : min(length_input, length_output))
                lhs(i:i) = rhs%characters(i)
            end do
            if (length_input < length_output) then
                do concurrent (i = length_input+1 : length_output)
                    lhs(i:i) = " "
                end do
            end if
        else
            do concurrent (i = 1 : length_output)
                lhs(i:i) = " "
            end do
        end if
    end subroutine

    elemental function concat_strings(lhs, rhs) result(concatenated)
        ! Sec. 3.3.2
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        type(varying_string) :: concatenated

        concatenated = char(lhs) // char(rhs)
    end function

    elemental function concat_string_and_character(lhs, rhs) result(concatenated)
        ! Sec. 3.3.2
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        type(varying_string) :: concatenated

        concatenated = char(lhs) // rhs
    end function

    elemental function concat_character_and_string(lhs, rhs) result(concatenated)
        ! Sec. 3.3.2
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        type(varying_string) :: concatenated

        concatenated = lhs // char(rhs)
    end function

    elemental function string_eq_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) == char(rhs)
    end function

    elemental function character_eq_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs == char(rhs)
    end function

    elemental function string_eq_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) == rhs
    end function

    elemental function string_ne_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) /= char(rhs)
    end function

    elemental function character_ne_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs /= char(rhs)
    end function

    elemental function string_ne_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) /= rhs
    end function

    elemental function string_lt_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) < char(rhs)
    end function

    elemental function character_lt_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs < char(rhs)
    end function

    elemental function string_lt_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) < rhs
    end function

    elemental function string_le_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) <= char(rhs)
    end function

    elemental function character_le_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs <= char(rhs)
    end function

    elemental function string_le_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) <= rhs
    end function

    elemental function string_gt_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) > char(rhs)
    end function

    elemental function character_gt_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs > char(rhs)
    end function

    elemental function string_gt_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) > rhs
    end function

    elemental function string_ge_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) >= char(rhs)
    end function

    elemental function character_ge_string(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        character(len=*), intent(in) :: lhs
        type(varying_string), intent(in) :: rhs
        logical :: equals

        equals = lhs >= char(rhs)
    end function

    elemental function string_ge_character(lhs, rhs) result(equals)
        ! Sec. 3.3.3
        type(varying_string), intent(in) :: lhs
        character(len=*), intent(in) :: rhs
        logical :: equals

        equals = char(lhs) >= rhs
    end function

    elemental function string_adjustl(string) result(adjusted)
        ! Sec. 3.4.1
        type(varying_string), intent(in) :: string
        type(varying_string) :: adjusted

        adjusted = adjustl(char(string))
    end function

    elemental function string_adjustr(string) result(adjusted)
        ! Sec. 3.4.2
        type(varying_string), intent(in) :: string
        type(varying_string) :: adjusted

        adjusted = adjustr(char(string))
    end function

    pure function string_to_char(string) result(chars)
        ! Sec. 3.4.3
        type(varying_string), intent(in) :: string
        character(len=:), allocatable :: chars

        if (allocated(string%characters)) then
            block
                character(len=size(string%characters)) :: tmp
                tmp = string
                chars = tmp
            end block
          else
            chars = ""
        end if
    end function

    pure function string_to_char_with_length(string, length) result(chars)
        ! Sec. 3.4.3
        type(varying_string), intent(in) :: string
        integer, intent(in) :: length
        character(len=max(0,length)) :: chars

        if (allocated(string%characters)) then
            chars = string
        end if
    end function

    elemental function string_iachar(c)
        ! Sec. 3.4.4
        type(varying_string), intent(in) :: c
        integer :: string_iachar

        string_iachar = iachar(char(c))
    end function

    elemental function string_ichar(c)
        ! Sec. 3.4.5
        type(varying_string), intent(in) :: c
        integer :: string_ichar

        string_ichar = ichar(char(c))
    end function

    elemental function string_index_string(string, substring, back) result(position)
        ! Sec. 3.4.6
        type(varying_string), intent(in) :: string
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: back
        integer :: position

        position = index(char(string), char(substring), back)
    end function

    elemental function string_index_character(string, substring, back) result(position)
        ! Sec. 3.4.6
        type(varying_string), intent(in) :: string
        character(len=*), intent(in) :: substring
        logical, optional, intent(in) :: back
        integer :: position

        position = index(char(string), substring, back)
    end function

    elemental function character_index_string(string, substring, back) result(position)
        ! Sec. 3.4.6
        character(len=*), intent(in) :: string
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: back
        integer :: position

        position = index(string, char(substring), back)
    end function

    elemental function len_string(string) result(length)
        ! Sec. 3.4.7
        type(varying_string), intent(in) :: string
        integer :: length

        length = len(char(string))
    end function

    elemental function len_trim_string(string) result(length)
        ! Sec. 3.4.8
        type(varying_string), intent(in) :: string
        integer :: length

        length = len_trim(char(string))
    end function

    elemental function string_lge_string(string_a, string_b) result(greater_than_or_equals)
        ! Sec 3.4.9
        type(varying_string), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: greater_than_or_equals

        greater_than_or_equals = lge(char(string_a), char(string_b))
    end function

    elemental function character_lge_string(string_a, string_b) result(greater_than_or_equals)
        ! Sec 3.4.9
        character(len=*), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: greater_than_or_equals

        greater_than_or_equals = lge(string_a, char(string_b))
    end function

    elemental function string_lge_character(string_a, string_b) result(greater_than_or_equals)
        ! Sec 3.4.9
        type(varying_string), intent(in) :: string_a
        character(len=*), intent(in) :: string_b
        logical :: greater_than_or_equals

        greater_than_or_equals = lge(char(string_a), string_b)
    end function

    elemental function string_lgt_string(string_a, string_b) result(greater_than)
        ! Sec 3.4.10
        type(varying_string), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: greater_than

        greater_than = lgt(char(string_a), char(string_b))
    end function

    elemental function character_lgt_string(string_a, string_b) result(greater_than)
        ! Sec 3.4.10
        character(len=*), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: greater_than

        greater_than = lgt(string_a, char(string_b))
    end function

    elemental function string_lgt_character(string_a, string_b) result(greater_than)
        ! Sec 3.4.10
        type(varying_string), intent(in) :: string_a
        character(len=*), intent(in) :: string_b
        logical :: greater_than

        greater_than = lgt(char(string_a), string_b)
    end function

    elemental function string_lle_string(string_a, string_b) result(less_than_or_equals)
        ! Sec 3.4.11
        type(varying_string), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: less_than_or_equals

        less_than_or_equals = lle(char(string_a), char(string_b))
    end function

    elemental function character_lle_string(string_a, string_b) result(less_than_or_equals)
        ! Sec 3.4.11
        character(len=*), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: less_than_or_equals

        less_than_or_equals = lle(string_a, char(string_b))
    end function

    elemental function string_lle_character(string_a, string_b) result(less_than_or_equals)
        ! Sec 3.4.11
        type(varying_string), intent(in) :: string_a
        character(len=*), intent(in) :: string_b
        logical :: less_than_or_equals

        less_than_or_equals = lle(char(string_a), string_b)
    end function

    elemental function string_llt_string(string_a, string_b) result(less_than)
        ! Sec 3.4.12
        type(varying_string), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        logical :: less_than
        intrinsic :: llt

        less_than = llt(char(string_a), char(string_b))
    end function

    elemental function character_llt_string(string_a, string_b) result(less_than)
        ! Sec 3.4.12
        character(len=*), intent(in) :: string_a
        type(varying_string), intent(in) :: string_b
        intrinsic :: llt
        logical :: less_than

        less_than = llt(string_a, char(string_b))
    end function

    elemental function string_llt_character(string_a, string_b) result(less_than)
        ! Sec 3.4.12
        type(varying_string), intent(in) :: string_a
        character(len=*), intent(in) :: string_b
        intrinsic :: llt
        logical :: less_than

        less_than = llt(char(string_a), string_b)
    end function
    elemental function string_repeat(string, ncopies) result(repeated)
        ! Sec. 3.4.13
        type(varying_string), intent(in) :: string
        integer, intent(in) :: ncopies
        type(varying_string) :: repeated

        intrinsic :: repeat

        if (allocated(string%characters)) then
            block
                character(len=len(string)) :: tmp_char
                tmp_char = string
                repeated = repeat(tmp_char, ncopies)
            end block
        else
            repeated = ""
        end if
    end function

    elemental function string_scan_string(string, set, back) result(position)
        ! Sec. 3.4.14
        type(varying_string), intent(in) :: string
        type(varying_string), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = scan(char(string), char(set), back)
    end function

    elemental function string_scan_character(string, set, back) result(position)
        ! Sec. 3.4.14
        type(varying_string), intent(in) :: string
        character(len=*), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = scan(char(string), set, back)
    end function

    elemental function character_scan_string(string, set, back) result(position)
        ! Sec. 3.4.14
        character(len=*), intent(in) :: string
        type(varying_string), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = scan(string, char(set), back)
    end function

    elemental function trim_string(string) result(trimmed)
        ! Sec. 3.4.15
        type(varying_string), intent(in) :: string
        type(varying_string) :: trimmed

        trimmed = trim(char(string))
    end function

    elemental function string_verify_string(string, set, back) result(position)
        ! Sec. 3.5.16
        type(varying_string), intent(in) :: string
        type(varying_string), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = verify(char(string), char(set), back)
    end function

    elemental function string_verify_character(string, set, back) result(position)
        ! Sec. 3.5.16
        type(varying_string), intent(in) :: string
        character(len=*), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = verify(char(string), set, back)
    end function

    elemental function character_verify_string(string, set, back) result(position)
        ! Sec. 3.5.16
        character(len=*), intent(in) :: string
        type(varying_string), intent(in) :: set
        logical, optional, intent(in) :: back
        integer :: position

        position = verify(string, char(set), back)
    end function

    elemental function var_str_char(char) result(var_str)
        ! Sec. 3.5.1
        character(len=*), intent(in) :: char
        type(varying_string) :: var_str

        var_str = char
    end function

    subroutine get_default_unit_to_end_of_record(string, maxlen, iostat)
        ! Sec. 3.6.1
        type(varying_string), intent(out) :: string
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        integer, parameter :: BUFFER_SIZE = 100
        character(len=BUFFER_SIZE) :: buffer
        integer :: next_read_length
        integer :: num_read
        integer :: num_to_read

        if (present(maxlen)) then
            num_to_read = maxlen
        else
            num_to_read = huge(1)
        end if
        string = ""
        if (present(iostat)) then
            do
                if (num_to_read <= 0) exit
                next_read_length = min(BUFFER_SIZE, num_to_read)
                read(*, fmt='(A)', advance='NO', eor=9999, size=num_read, iostat=iostat) buffer(1:next_read_length)
                if (iostat /= 0) return
                string = string // buffer(1:next_read_length)
                num_to_read = num_to_read - next_read_length
            end do
        else
            do
                if (num_to_read <= 0) exit
                next_read_length = min(BUFFER_SIZE, num_to_read)
                read(*, fmt='(A)', advance='NO', eor=9999, size=num_read) buffer(1:next_read_length)
                string = string // buffer(1:next_read_length)
                num_to_read = num_to_read - next_read_length
            end do
        end if
        return
        9999 string = string // buffer(1:num_read)
    end subroutine

    subroutine get_with_unit_to_end_of_record(unit, string, maxlen, iostat)
        ! Sec. 3.6.1
        integer, intent(in) :: unit
        type(varying_string), intent(out) :: string
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        integer, parameter :: BUFFER_SIZE = 100
        character(len=BUFFER_SIZE) :: buffer
        integer :: next_read_length
        integer :: num_read
        integer :: num_to_read

        if (present(maxlen)) then
            num_to_read = maxlen
        else
            num_to_read = huge(1)
        end if
        string = ""
        if (present(iostat)) then
            do
                if (num_to_read <= 0) exit
                next_read_length = min(BUFFER_SIZE, num_to_read)
                read(unit, fmt='(A)', advance='NO', eor=9999, size=num_read, iostat=iostat) buffer(1:next_read_length)
                if (iostat /= 0) return
                string = string // buffer(1:next_read_length)
                num_to_read = num_to_read - next_read_length
            end do
        else
            do
                if (num_to_read <= 0) exit
                next_read_length = min(BUFFER_SIZE, num_to_read)
                read(unit, fmt='(A)', advance='NO', eor=9999, size=num_read) buffer(1:next_read_length)
                string = string // buffer(1:next_read_length)
                num_to_read = num_to_read - next_read_length
            end do
        end if
        return
        9999 string = string // buffer(1:num_read)
    end subroutine

    subroutine get_default_unit_to_terminator_string(string, set, separator, maxlen, iostat)
        ! Sec. 3.6.1
        type(varying_string), intent(out) :: string
        type(varying_string), intent(in) :: set ! possible terminator characters
        type(varying_string), optional, intent(out) :: separator ! actual terminator
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        call get(string, char(set), separator, maxlen, iostat)
    end subroutine

    subroutine get_with_unit_to_terminator_string(unit, string, set, separator, maxlen, iostat)
        ! Sec. 3.6.1
        integer, intent(in) :: unit
        type(varying_string), intent(out) :: string
        type(varying_string), intent(in) :: set ! possible terminator characters
        type(varying_string), optional, intent(out) :: separator ! actual terminator
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        call get(unit, string, char(set), separator, maxlen, iostat)
    end subroutine

    subroutine get_default_unit_to_terminator_characters(string, set, separator, maxlen, iostat)
        ! Sec. 3.6.1
        type(varying_string), intent(out) :: string
        character(len=*), intent(in) :: set ! possible terminator characters
        type(varying_string), optional, intent(out) :: separator ! actual terminator
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        character(len=1) :: buffer
        integer :: num_to_read

        if (present(maxlen)) then
            num_to_read = maxlen
        else
            num_to_read = huge(1)
        end if
        string = ""
        if (present(separator)) separator = ""
        if (present(iostat)) then
            do
                if (num_to_read <= 0) exit
                read(*, fmt='(A)', advance='NO', eor=9999, iostat=iostat) buffer
                if (iostat /= 0) return
                if (index(set, buffer) /= 0) then
                    if (present(separator)) separator = buffer
                    return
                end if
                string = string // buffer
                num_to_read = num_to_read - 1
            end do
        else
            do
                if (num_to_read <= 0) exit
                read(*, fmt='(A)', advance='NO', eor=9999) buffer
                if (index(set, buffer) /= 0) then
                    if (present(separator)) separator = buffer
                    return
                end if
                string = string // buffer
                num_to_read = num_to_read - 1
            end do
        end if
        9999 continue
    end subroutine

    subroutine get_with_unit_to_terminator_characters(unit, string, set, separator, maxlen, iostat)
        ! Sec. 3.6.1
        integer, intent(in) :: unit
        type(varying_string), intent(out) :: string
        character(len=*), intent(in) :: set ! possible terminator characters
        type(varying_string), optional, intent(out) :: separator ! actual terminator
        integer, optional, intent(in) :: maxlen
        integer, optional, intent(out) :: iostat

        character(len=1) :: buffer
        integer :: num_to_read

        if (present(maxlen)) then
            num_to_read = maxlen
        else
            num_to_read = huge(1)
        end if
        string = ""
        if (present(separator)) separator = ""
        if (present(iostat)) then
            do
                if (num_to_read <= 0) exit
                read(unit, fmt='(A)', advance='NO', eor=9999, iostat=iostat) buffer
                if (iostat /= 0) return
                if (index(set, buffer) /= 0) then
                    if (present(separator)) separator = buffer
                    return
                end if
                string = string // buffer
                num_to_read = num_to_read - 1
            end do
        else
            do
                if (num_to_read <= 0) exit
                read(unit, fmt='(A)', advance='NO', eor=9999) buffer
                if (index(set, buffer) /= 0) then
                    if (present(separator)) separator = buffer
                    return
                end if
                string = string // buffer
                num_to_read = num_to_read - 1
            end do
        end if
        9999 continue
    end subroutine

    subroutine put_string_default_unit(string, iostat)
        ! Sec. 3.6.2
        type(varying_string), intent(in) :: string
        integer, optional, intent(out) :: iostat

        call put(char(string), iostat)
    end subroutine

    subroutine put_string_with_unit(unit, string, iostat)
        ! Sec. 3.6.2
        integer, intent(in) :: unit
        type(varying_string), intent(in) :: string
        integer, optional, intent(out) :: iostat

        call put(unit, char(string), iostat)
    end subroutine

    subroutine put_characters_default_unit(string, iostat)
        ! Sec. 3.6.2
        character(len=*), intent(in) :: string
        integer, optional, intent(out) :: iostat

        if (present(iostat)) then
            write(*, fmt='(A)', advance='NO', iostat=iostat) string
        else
            write(*, fmt='(A)', advance='NO') string
        end if
    end subroutine

    subroutine put_characters_with_unit(unit, string, iostat)
        ! Sec. 3.6.2
        integer, intent(in) :: unit
        character(len=*), intent(in) :: string
        integer, optional, intent(out) :: iostat

        if (present(iostat)) then
            write(unit, fmt='(A)', advance='NO', iostat=iostat) string
        else
            write(unit, fmt='(A)', advance='NO') string
        end if
    end subroutine

    subroutine put_line_string_default_unit(string, iostat)
        ! Sec. 3.6.3
        type(varying_string), intent(in) :: string
        integer, optional, intent(out) :: iostat

        call put_line(char(string), iostat)
    end subroutine

    subroutine put_line_string_with_unit(unit, string, iostat)
        ! Sec. 3.6.3
        integer, intent(in) :: unit
        type(varying_string), intent(in) :: string
        integer, optional, intent(out) :: iostat

        call put_line(unit, char(string), iostat)
    end subroutine

    subroutine put_line_characters_default_unit(string, iostat)
        ! Sec. 3.6.3
        character(len=*), intent(in) :: string
        integer, optional, intent(out) :: iostat

        if (present(iostat)) then
            write(*, fmt='(A,/)', advance='NO', iostat=iostat) string
        else
            write(*, fmt='(A,/)', advance='NO') string
        end if
    end subroutine

    subroutine put_line_characters_with_unit(unit, string, iostat)
        ! Sec. 3.6.3
        integer, intent(in) :: unit
        character(len=*), intent(in) :: string
        integer, optional, intent(out) :: iostat

        if (present(iostat)) then
            write(unit, fmt='(A,/)', advance='NO', iostat=iostat) string
        else
            write(unit, fmt='(A,/)', advance='NO') string
        end if
    end subroutine

    elemental function extract_character(string, start, finish) result(extracted)
        ! Sec. 3.7.1
        character(len=*), intent(in) :: string
        integer, optional, intent(in) :: start
        integer, optional, intent(in) :: finish
        type(varying_string) :: extracted

        integer :: start_
        integer :: finish_

        if (present(start)) then
            start_ = max(1, start)
        else
            start_ = 1
        end if
        if (present(finish)) then
            finish_ = min(len(string), finish)
        else
            finish_ = len(string)
        end if

        extracted = string(start_:finish_)
    end function

    elemental function extract_string(string, start, finish) result(extracted)
        ! Sec. 3.7.1
        type(varying_string), intent(in) :: string
        integer, optional, intent(in) :: start
        integer, optional, intent(in) :: finish
        type(varying_string) :: extracted

        extracted = extract(char(string), start, finish)
    end function

    elemental function insert_character_into_character(string, start, substring) result(inserted)
        ! Sec. 3.7.2
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        character(len=*), intent(in) :: substring
        type(varying_string) :: inserted

        type(varying_string) :: beginning
        type(varying_string) :: middle
        type(varying_string) :: end_

        if (start <= 1) then
            beginning = substring
            middle = string
            end_ = ""
        else if (start > len(string)) then
            beginning = string
            middle = substring
            end_ = ""
        else
            beginning = string(1:start-1)
            middle = substring
            end_ = string(start:)
        end if
        inserted = beginning // middle // end_
    end function

    elemental function insert_character_into_string(string, start, substring) result(inserted)
        ! Sec. 3.7.2
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        character(len=*), intent(in) :: substring
        type(varying_string) :: inserted

        inserted = insert(char(string), start, substring)
    end function

    elemental function insert_string_into_character(string ,start, substring) result(inserted)
        ! Sec. 3.7.2
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        type(varying_string), intent(in) :: substring
        type(varying_string) :: inserted

        inserted = insert(string, start, char(substring))
    end function

    elemental function insert_string_into_string(string ,start, substring) result(inserted)
        ! Sec. 3.7.2
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        type(varying_string), intent(in) :: substring
        type(varying_string) :: inserted

        inserted = insert(char(string), start, char(substring))
    end function

    elemental function remove_character(string, start, finish) result(removed)
        ! Sec. 3.7.3
        character(len=*), intent(in) :: string
        integer, optional, intent(in) :: start
        integer, optional, intent(in) :: finish
        type(varying_string) :: removed

        integer :: start_
        integer :: finish_
        type(varying_string) :: beginning
        type(varying_string) :: end_

        if (present(start)) then
            start_ = start
        else
            start_ = 1
        end if
        if (present(finish)) then
            finish_ = finish
        else
            finish_ = len(string)
        end if

        if (start_ > finish_) then
            removed = string
        else
            beginning = string(1:start_ - 1)
            end_ = string(finish_ + 1:len(string))
            removed = beginning // end_
        end if
    end function

    elemental function remove_string(string, start, finish) result(removed)
        ! Sec. 3.7.3
        type(varying_string), intent(in) :: string
        integer, optional, intent(in) :: start
        integer, optional, intent(in) :: finish
        type(varying_string) :: removed

        removed = remove(char(string), start, finish)
    end function

    elemental function replace_character_with_character_start( &
                string, start, substring) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        character(len=*), intent(in) :: substring
        type(varying_string) :: replaced

        integer :: start_

        start_ = max(1, start)
        replaced = insert( &
                remove(string, start_, start_ + len(substring) - 1), &
                start_, &
                substring)
    end function

    elemental function replace_string_with_character_start( &
            string, start, substring) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        character(len=*), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(char(string), start, substring)
    end function

    elemental function replace_character_with_string_start( &
            string, start, substring) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        type(varying_string), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(string, start, char(substring))
    end function

    elemental function replace_string_with_string_start( &
            string, start, substring) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        type(varying_string), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(char(string), start, char(substring))
    end function

    elemental function replace_character_with_character_range( &
            string, start, finish, substring) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        integer, intent(in) :: finish
        character(len=*), intent(in) :: substring
        type(varying_string) :: replaced

        type(varying_string) :: beginning
        type(varying_string) :: ending

        beginning = string(1 : start-1)
        ending = string(max(finish+1, start) : )
        replaced = beginning // substring // ending
    end function

    elemental function replace_string_with_character_range( &
            string, start, finish, substring) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        integer, intent(in) :: finish
        character(len=*), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(char(string), start, finish, substring)
    end function

    elemental function replace_character_with_string_range( &
            string, start, finish, substring) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        integer, intent(in) :: start
        integer, intent(in) :: finish
        type(varying_string), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(string, start, finish, char(substring))
    end function

    elemental function replace_string_with_string_range( &
            string, start, finish, substring) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        integer, intent(in) :: start
        integer, intent(in) :: finish
        type(varying_string), intent(in) :: substring
        type(varying_string) :: replaced

        replaced = replace(char(string), start, finish, char(substring))
    end function

    elemental function replace_target_character_with_character_in_character( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        character(len=*), intent(in) :: target
        character(len=*), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        logical :: back_
        logical :: every_

        if (present(back)) then
            back_ = back
        else
            back_ = .false.
        end if

        if (present(every)) then
            every_ = every
        else
            every_ = .false.
        end if

        if (every_) then
            block
                integer :: i
                integer :: new_length
                character(len=:), allocatable :: new_string
                integer :: num_targets
                integer :: target_length, substring_length
                integer, allocatable :: target_positions(:), substring_positions(:)

                num_targets = get_num_targets(string, target, back_)
                if (num_targets > 0) then
                    target_length = len(target)
                    substring_length = len(substring)
                    new_length = &
                            len(string) &
                            - num_targets*target_length &
                            + num_targets*substring_length
                    allocate(character(len=new_length)::new_string)
                    allocate(target_positions(num_targets))
                    allocate(substring_positions(num_targets))
                    call get_positions( &
                            string, &
                            target, &
                            back_, &
                            num_targets, &
                            substring_length, &
                            target_positions, &
                            substring_positions)
                    if (target_positions(1) > 1) then
                        new_string(1 : substring_positions(1)-1) = string(1 : target_positions(1)-1)
                    end if
                    new_string(substring_positions(1) : substring_positions(1)+substring_length-1) = substring
                    do i = 2, num_targets
                        new_string(substring_positions(i-1)+substring_length : substring_positions(i)-1) = &
                                string(target_positions(i-1)+target_length:target_positions(i)-1)
                        new_string(substring_positions(i) : substring_positions(i)+substring_length-1) = substring
                    end do
                    if (target_positions(num_targets) + target_length <= len(string)) then
                        new_string(substring_positions(num_targets)+substring_length:) = &
                                string(target_positions(num_targets)+target_length:)
                    end if
                    replaced = new_string
                else
                    replaced = string
                end if
            end block
        else
            block
                integer :: position

                position = index(string, target, back_)
                if (position /= 0) then
                    replaced = string(1:position-1) // substring // string(position+len(target):)
                else
                    replaced = string
                end if
            end block
        end if
    contains
        pure function get_num_targets(string_, target_, back__) result(num_targets_)
            character(len=*), intent(in) :: string_
            character(len=*), intent(in) :: target_
            logical, intent(in) :: back__
            integer :: num_targets_

            integer :: len_string
            integer :: len_target
            integer :: prev_start

            len_string = len(string_)
            len_target = len(target_)
            prev_start = index(string_, target_, back__)
            if (prev_start == 0) then
                num_targets_ = 0
            else
                num_targets_ = 1
                if (back__) then
                  do
                      if (prev_start == 1) exit
                      prev_start = index(string_(1:prev_start-1), target_, back__)
                      if (prev_start == 0) then
                          exit
                      else
                          num_targets_ = num_targets_ + 1
                      end if
                  end do
                else
                    do
                        if (prev_start+len_target > len_string) exit
                        if (index(string_(prev_start+len_target:), target_) == 0) then
                            exit
                        else
                            prev_start = index(string_(prev_start+len_target:), target_) + prev_start + len_target - 1
                            num_targets_ = num_targets_ + 1
                        end if
                    end do
                end if
            end if
        end function

        pure subroutine get_positions( &
              string_, &
              target_, &
              back__, &
              num_targets_, &
              len_substring, &
              target_positions_, &
              substring_positions_)
          character(len=*), intent(in) :: string_
          character(len=*), intent(in) :: target_
          logical, intent(in) :: back__
          integer, intent(in) :: num_targets_
          integer, intent(in) :: len_substring
          integer, intent(out) :: target_positions_(:), substring_positions_(:)

          integer :: i, len_target

          len_target = len(target_)
          if (back__) then
              target_positions_(num_targets_) = index(string_, target_, back__)
              do i = num_targets_-1, 1, -1
                  target_positions_(i) = index(string_(1:target_positions_(i+1)-1), target_, back__)
              end do
          else
              target_positions_(1) = index(string_, target_)
              do i = 2, num_targets_
                  target_positions_(i) = &
                          index(string_(target_positions_(i-1)+len_target:), target_) &
                          + target_positions_(i-1) + len_target - 1
              end do
          end if
          substring_positions_(1) = target_positions_(1)
          do i = 2, num_targets_
              substring_positions_(i) = &
                      substring_positions_(i-1) &
                      + (target_positions_(i) - target_positions_(i-1)) &
                      + (len_substring - len_target)
          end do
        end subroutine
    end function replace_target_character_with_character_in_character

    elemental function replace_target_character_with_character_in_string( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        character(len=*), intent(in) :: target
        character(len=*), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(char(string), target, substring, every, back)
    end function

    elemental function replace_target_character_with_string_in_character( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        character(len=*), intent(in) :: target
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(string, target, char(substring), every, back)
    end function

    elemental function replace_target_character_with_string_in_string( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        character(len=*), intent(in) :: target
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(char(string), target, char(substring), every, back)
    end function

    elemental function replace_target_string_with_character_in_character( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        type(varying_string), intent(in) :: target
        character(len=*), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(string, char(target), substring, every, back)
    end function

    elemental function replace_target_string_with_character_in_string( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        type(varying_string), intent(in) :: target
        character(len=*), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(char(string), char(target), substring, every, back)
    end function

    elemental function replace_target_string_with_string_in_character( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        character(len=*), intent(in) :: string
        type(varying_string), intent(in) :: target
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(string, char(target), char(substring), every, back)
    end function

    elemental function replace_target_string_with_string_in_string( &
            string, target, substring, every, back) result(replaced)
        ! Sec. 3.7.4
        type(varying_string), intent(in) :: string
        type(varying_string), intent(in) :: target
        type(varying_string), intent(in) :: substring
        logical, optional, intent(in) :: every
        logical, optional, intent(in) :: back
        type(varying_string) :: replaced

        replaced = replace(char(string), char(target), char(substring), every, back)
    end function

    elemental subroutine split_character(string, word, set, separator, back)
        ! Sec. 3.7.5
        type(varying_string), intent(inout) :: string
        type(varying_string), intent(out) :: word
        character(len=*), intent(in) :: set
        type(varying_string), optional, intent(out) :: separator
        logical, optional, intent(in) :: back

        logical :: backwards
        integer :: i
        integer :: string_length
        character(len=:), allocatable :: temp

        allocate(character(len=0) :: temp) ! TODO: remove once gfortran bug is fixed
        string_length = len(string)
        if (present(back)) then
            backwards = back
        else
            backwards = .false.
        end if
        if (backwards) then
            do i = string_length, 1, -1
                if (index(set, extract(string, i, i)) /= 0) exit
            end do
            if (i < 1) then
                word = string
                string = ""
                if (present(separator)) separator = ""
            else
                word = extract(string, i+1)
                temp = char(extract(string, 1, i-1))
                if (present(separator)) separator = extract(string, i, i)
                string = temp
            end if
        else
            do i = 1, string_length
                if (index(set, extract(string, i, i)) /= 0) exit
            end do
            if (i > string_length) then
                word = string
                string = ""
                if (present(separator)) separator = ""
            else
                word = extract(string, 1, i-1)
                temp = char(extract(string, i+1))
                if (present(separator)) separator = extract(string, i, i)
                string = temp
            end if
        end if
    end subroutine

    elemental subroutine split_string(string, word, set, separator, back)
        ! Sec. 3.7.5
        type(varying_string), intent(inout) :: string
        type(varying_string), intent(out) :: word
        type(varying_string), intent(in) :: set
        type(varying_string), optional, intent(out) :: separator
        logical, optional, intent(in) :: back

        call split(string, word, char(set), separator, back)
    end subroutine

    subroutine write_formatted(self, unit, iotype, v_list, iostat, iomsg)
        class(varying_string), intent(in) :: self
        integer, intent(in) :: unit
        character(len=*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(len=*), intent(inout) :: iomsg

        select case (iotype)
        case ("LISTDIRECTED")
            write(unit, '(A)', iostat=iostat, iomsg=iomsg) char(self)
        case ("NAMELIST")
            write(unit, '(A,",")', iostat=iostat, iomsg=iomsg) '"' // char(self) // '"'
        case ("DT")
            if (size(v_list) > 0) then
                if (v_list(1) > 0) then
                    block
                        character(len=20) :: num_string
                        character(len=:), allocatable :: fmt

                        write(num_string, '(I0)') v_list(1)
                        fmt = "(A" // trim(num_string) // ")"
                        write(unit, fmt, iostat=iostat, iomsg=iomsg) char(self)
                    end block
                else
                    write(unit, '(A)', iostat=iostat, iomsg=iomsg) char(self)
                end if
            else
                write(unit, '(A)', iostat=iostat, iomsg=iomsg) char(self)
            end if
        case default
            error stop "Unknown iotype: " // iotype
        end select
    end subroutine

    subroutine write_unformatted(self, unit, iostat, iomsg)
        class(varying_string), intent(in) :: self
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len=*), intent(inout) :: iomsg

        write(unit, iostat=iostat, iomsg=iomsg) char(self)
    end subroutine

    subroutine read_formatted(self, unit, iotype, v_list, iostat, iomsg)
        class(varying_string), intent(inout) :: self
        integer, intent(in) :: unit
        character(len=*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(len=*), intent(inout) :: iomsg

        select type (self)
        type is (varying_string)
        select case (iotype)
        case ("LISTDIRECTED", "NAMELIST")
            block
                character(len=1) :: next_char
                character(len=:), allocatable :: end_chars
                next_char = " "
                do while (index(" ", next_char) > 0)
                    read(unit, fmt='(A)', iostat=iostat, iomsg=iomsg) next_char
                    if (iostat /= 0) return
                end do
                if (next_char == ',') then
                    return
                else if (next_char == '"') then
                    end_chars = '"'
                    self = ""
                else if (next_char == "'") then
                    end_chars = "'"
                    self = ""
                else
                    end_chars = " ,"
                    self = next_char
                end if
                do
                    read(unit, fmt='(A)', iostat=iostat, iomsg=iomsg) next_char
                    if (iostat /= 0) return
                    if (index(end_chars, next_char) > 0) return
                    self = self // next_char
                end do
            end block
        case ("DT")
            if (size(v_list) > 0) then
                if (v_list(1) > 0) then
                    block
                        character(len=v_list(1)) :: to_read
                        read(unit, '(A)', iostat=iostat, iomsg=iomsg) to_read
                        self = to_read
                    end block
                else
                    read(unit, '(DT)', iostat=iostat, iomsg=iomsg) self
                end if
            else
                block
                    character(len=100) :: buffer
                    integer :: num_read
                    logical :: first
                    first = .true.
                    do
                        read(unit, fmt='(A)', advance='no', iostat=iostat, iomsg=iomsg, size=num_read) buffer
                        if (is_iostat_eor(iostat) .or. is_iostat_end(iostat)) then
                            if (first) then
                                self = buffer(1:num_read)
                            else
                                self = self // buffer(1:num_read)
                            end if
                            return
                        else if (iostat /= 0) then
                            return
                        else
                            if (first) then
                                self = buffer
                                first = .false.
                            else
                                self = self // buffer
                            end if
                        end if
                    end do
                end block
            end if
        case default
            error stop "Unknown iotype: " // iotype
        end select
        class default
            error stop "Not supported for extended types"
        end select
    end subroutine

    subroutine read_unformatted(self, unit, iostat, iomsg)
        class(varying_string), intent(inout) :: self
        integer, intent(in) :: unit
        integer, intent(out) :: iostat
        character(len=*), intent(inout) :: iomsg

        read(unit, *, iostat=iostat, iomsg=iomsg) self
    end subroutine
end module
