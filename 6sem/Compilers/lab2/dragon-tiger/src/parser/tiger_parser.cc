// A Bison parser, made by GNU Bison 3.5.1.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2020 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// Undocumented macros, especially those whose name start with YY_,
// are private implementation details.  Do not rely on them.





#include "tiger_parser.hh"


// Unqualified %code blocks.
#line 35 "tiger_parser.yy"

#include "parser_driver.hh"

#line 49 "tiger_parser.cc"


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

namespace yy {
#line 140 "tiger_parser.cc"


  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  tiger_parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr;
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              else
                goto append;

            append:
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }


  /// Build a parser object.
  tiger_parser::tiger_parser (ParserDriver& driver_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      driver (driver_yyarg)
  {}

  tiger_parser::~tiger_parser ()
  {}

  tiger_parser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------------.
  | Symbol types.  |
  `---------------*/



  // by_state.
  tiger_parser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  tiger_parser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  tiger_parser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  tiger_parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  tiger_parser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  tiger_parser::symbol_number_type
  tiger_parser::by_state::type_get () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[+state];
  }

  tiger_parser::stack_symbol_type::stack_symbol_type ()
  {}

  tiger_parser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.location))
  {
    switch (that.type_get ())
    {
      case 44: // decl
      case 46: // varDecl
      case 47: // funcDecl
        value.YY_MOVE_OR_COPY< Decl * > (YY_MOVE (that.value));
        break;

      case 43: // program
      case 45: // expr
      case 48: // stringExpr
      case 49: // intExpr
      case 50: // var
      case 51: // callExpr
      case 52: // ifExpr
      case 53: // negExpr
      case 54: // opExpr
      case 55: // assignExpr
      case 56: // whileExpr
      case 57: // forExpr
      case 58: // breakExpr
      case 59: // letExpr
      case 60: // seqExpr
        value.YY_MOVE_OR_COPY< Expr * > (YY_MOVE (that.value));
        break;

      case 37: // "id"
      case 38: // "string"
        value.YY_MOVE_OR_COPY< Symbol > (YY_MOVE (that.value));
        break;

      case 68: // param
        value.YY_MOVE_OR_COPY< VarDecl * > (YY_MOVE (that.value));
        break;

      case 69: // typeannotation
        value.YY_MOVE_OR_COPY< boost::optional<Symbol> > (YY_MOVE (that.value));
        break;

      case 39: // "int"
        value.YY_MOVE_OR_COPY< int > (YY_MOVE (that.value));
        break;

      case 67: // decls
        value.YY_MOVE_OR_COPY< std::vector<Decl *> > (YY_MOVE (that.value));
        break;

      case 61: // exprs
      case 62: // nonemptyexprs
      case 63: // arguments
      case 64: // nonemptyarguments
        value.YY_MOVE_OR_COPY< std::vector<Expr *> > (YY_MOVE (that.value));
        break;

      case 65: // params
      case 66: // nonemptyparams
        value.YY_MOVE_OR_COPY< std::vector<VarDecl *> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  tiger_parser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.location))
  {
    switch (that.type_get ())
    {
      case 44: // decl
      case 46: // varDecl
      case 47: // funcDecl
        value.move< Decl * > (YY_MOVE (that.value));
        break;

      case 43: // program
      case 45: // expr
      case 48: // stringExpr
      case 49: // intExpr
      case 50: // var
      case 51: // callExpr
      case 52: // ifExpr
      case 53: // negExpr
      case 54: // opExpr
      case 55: // assignExpr
      case 56: // whileExpr
      case 57: // forExpr
      case 58: // breakExpr
      case 59: // letExpr
      case 60: // seqExpr
        value.move< Expr * > (YY_MOVE (that.value));
        break;

      case 37: // "id"
      case 38: // "string"
        value.move< Symbol > (YY_MOVE (that.value));
        break;

      case 68: // param
        value.move< VarDecl * > (YY_MOVE (that.value));
        break;

      case 69: // typeannotation
        value.move< boost::optional<Symbol> > (YY_MOVE (that.value));
        break;

      case 39: // "int"
        value.move< int > (YY_MOVE (that.value));
        break;

      case 67: // decls
        value.move< std::vector<Decl *> > (YY_MOVE (that.value));
        break;

      case 61: // exprs
      case 62: // nonemptyexprs
      case 63: // arguments
      case 64: // nonemptyarguments
        value.move< std::vector<Expr *> > (YY_MOVE (that.value));
        break;

      case 65: // params
      case 66: // nonemptyparams
        value.move< std::vector<VarDecl *> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

    // that is emptied.
    that.type = empty_symbol;
  }

#if YY_CPLUSPLUS < 201103L
  tiger_parser::stack_symbol_type&
  tiger_parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    switch (that.type_get ())
    {
      case 44: // decl
      case 46: // varDecl
      case 47: // funcDecl
        value.copy< Decl * > (that.value);
        break;

      case 43: // program
      case 45: // expr
      case 48: // stringExpr
      case 49: // intExpr
      case 50: // var
      case 51: // callExpr
      case 52: // ifExpr
      case 53: // negExpr
      case 54: // opExpr
      case 55: // assignExpr
      case 56: // whileExpr
      case 57: // forExpr
      case 58: // breakExpr
      case 59: // letExpr
      case 60: // seqExpr
        value.copy< Expr * > (that.value);
        break;

      case 37: // "id"
      case 38: // "string"
        value.copy< Symbol > (that.value);
        break;

      case 68: // param
        value.copy< VarDecl * > (that.value);
        break;

      case 69: // typeannotation
        value.copy< boost::optional<Symbol> > (that.value);
        break;

      case 39: // "int"
        value.copy< int > (that.value);
        break;

      case 67: // decls
        value.copy< std::vector<Decl *> > (that.value);
        break;

      case 61: // exprs
      case 62: // nonemptyexprs
      case 63: // arguments
      case 64: // nonemptyarguments
        value.copy< std::vector<Expr *> > (that.value);
        break;

      case 65: // params
      case 66: // nonemptyparams
        value.copy< std::vector<VarDecl *> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }

  tiger_parser::stack_symbol_type&
  tiger_parser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    switch (that.type_get ())
    {
      case 44: // decl
      case 46: // varDecl
      case 47: // funcDecl
        value.move< Decl * > (that.value);
        break;

      case 43: // program
      case 45: // expr
      case 48: // stringExpr
      case 49: // intExpr
      case 50: // var
      case 51: // callExpr
      case 52: // ifExpr
      case 53: // negExpr
      case 54: // opExpr
      case 55: // assignExpr
      case 56: // whileExpr
      case 57: // forExpr
      case 58: // breakExpr
      case 59: // letExpr
      case 60: // seqExpr
        value.move< Expr * > (that.value);
        break;

      case 37: // "id"
      case 38: // "string"
        value.move< Symbol > (that.value);
        break;

      case 68: // param
        value.move< VarDecl * > (that.value);
        break;

      case 69: // typeannotation
        value.move< boost::optional<Symbol> > (that.value);
        break;

      case 39: // "int"
        value.move< int > (that.value);
        break;

      case 67: // decls
        value.move< std::vector<Decl *> > (that.value);
        break;

      case 61: // exprs
      case 62: // nonemptyexprs
      case 63: // arguments
      case 64: // nonemptyarguments
        value.move< std::vector<Expr *> > (that.value);
        break;

      case 65: // params
      case 66: // nonemptyparams
        value.move< std::vector<VarDecl *> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  tiger_parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  tiger_parser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
#if defined __GNUC__ && ! defined __clang__ && ! defined __ICC && __GNUC__ * 100 + __GNUC_MINOR__ <= 408
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
#endif
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    YYUSE (yytype);
    yyo << ')';
  }
#endif

  void
  tiger_parser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  tiger_parser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  tiger_parser::yypop_ (int n)
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  tiger_parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  tiger_parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  tiger_parser::debug_level_type
  tiger_parser::debug_level () const
  {
    return yydebug_;
  }

  void
  tiger_parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  tiger_parser::state_type
  tiger_parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  bool
  tiger_parser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  tiger_parser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  tiger_parser::operator() ()
  {
    return parse ();
  }

  int
  tiger_parser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    // User initialization code.
#line 26 "tiger_parser.yy"
{
  // Initialize the initial location.
  yyla.location.begin.filename = yyla.location.end.filename = &driver.file;
}

#line 678 "tiger_parser.cc"


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            symbol_type yylookahead (yylex (driver));
            yyla.move (yylookahead);
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      {
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
      switch (yyr1_[yyn])
    {
      case 44: // decl
      case 46: // varDecl
      case 47: // funcDecl
        yylhs.value.emplace< Decl * > ();
        break;

      case 43: // program
      case 45: // expr
      case 48: // stringExpr
      case 49: // intExpr
      case 50: // var
      case 51: // callExpr
      case 52: // ifExpr
      case 53: // negExpr
      case 54: // opExpr
      case 55: // assignExpr
      case 56: // whileExpr
      case 57: // forExpr
      case 58: // breakExpr
      case 59: // letExpr
      case 60: // seqExpr
        yylhs.value.emplace< Expr * > ();
        break;

      case 37: // "id"
      case 38: // "string"
        yylhs.value.emplace< Symbol > ();
        break;

      case 68: // param
        yylhs.value.emplace< VarDecl * > ();
        break;

      case 69: // typeannotation
        yylhs.value.emplace< boost::optional<Symbol> > ();
        break;

      case 39: // "int"
        yylhs.value.emplace< int > ();
        break;

      case 67: // decls
        yylhs.value.emplace< std::vector<Decl *> > ();
        break;

      case 61: // exprs
      case 62: // nonemptyexprs
      case 63: // arguments
      case 64: // nonemptyarguments
        yylhs.value.emplace< std::vector<Expr *> > ();
        break;

      case 65: // params
      case 66: // nonemptyparams
        yylhs.value.emplace< std::vector<VarDecl *> > ();
        break;

      default:
        break;
    }


      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 2:
#line 122 "tiger_parser.yy"
              { driver.result_ast = yystack_[0].value.as < Expr * > (); }
#line 862 "tiger_parser.cc"
    break;

  case 3:
#line 125 "tiger_parser.yy"
              { yylhs.value.as < Decl * > () = yystack_[0].value.as < Decl * > (); }
#line 868 "tiger_parser.cc"
    break;

  case 4:
#line 126 "tiger_parser.yy"
              { yylhs.value.as < Decl * > () = yystack_[0].value.as < Decl * > (); }
#line 874 "tiger_parser.cc"
    break;

  case 5:
#line 129 "tiger_parser.yy"
                 { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 880 "tiger_parser.cc"
    break;

  case 6:
#line 130 "tiger_parser.yy"
             { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 886 "tiger_parser.cc"
    break;

  case 7:
#line 131 "tiger_parser.yy"
             { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 892 "tiger_parser.cc"
    break;

  case 8:
#line 132 "tiger_parser.yy"
         { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 898 "tiger_parser.cc"
    break;

  case 9:
#line 133 "tiger_parser.yy"
              { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 904 "tiger_parser.cc"
    break;

  case 10:
#line 134 "tiger_parser.yy"
            { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 910 "tiger_parser.cc"
    break;

  case 11:
#line 135 "tiger_parser.yy"
             { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 916 "tiger_parser.cc"
    break;

  case 12:
#line 136 "tiger_parser.yy"
            { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 922 "tiger_parser.cc"
    break;

  case 13:
#line 137 "tiger_parser.yy"
                { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 928 "tiger_parser.cc"
    break;

  case 14:
#line 138 "tiger_parser.yy"
               { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 934 "tiger_parser.cc"
    break;

  case 15:
#line 139 "tiger_parser.yy"
             { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 940 "tiger_parser.cc"
    break;

  case 16:
#line 140 "tiger_parser.yy"
               { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 946 "tiger_parser.cc"
    break;

  case 17:
#line 141 "tiger_parser.yy"
             { yylhs.value.as < Expr * > () = yystack_[0].value.as < Expr * > (); }
#line 952 "tiger_parser.cc"
    break;

  case 18:
#line 145 "tiger_parser.yy"
  { yylhs.value.as < Decl * > () = new VarDecl(yystack_[4].location, yystack_[3].value.as < Symbol > (), yystack_[2].value.as < boost::optional<Symbol> > (), yystack_[0].value.as < Expr * > ()); }
#line 958 "tiger_parser.cc"
    break;

  case 19:
#line 149 "tiger_parser.yy"
  { yylhs.value.as < Decl * > () = new FunDecl(yystack_[7].location, yystack_[6].value.as < Symbol > (), yystack_[2].value.as < boost::optional<Symbol> > (), yystack_[4].value.as < std::vector<VarDecl *> > (), yystack_[0].value.as < Expr * > ()); }
#line 964 "tiger_parser.cc"
    break;

  case 20:
#line 155 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new StringLiteral(yystack_[0].location, yystack_[0].value.as < Symbol > ()); }
#line 970 "tiger_parser.cc"
    break;

  case 21:
#line 159 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new IntegerLiteral(yystack_[0].location, yystack_[0].value.as < int > ()); }
#line 976 "tiger_parser.cc"
    break;

  case 22:
#line 163 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new Identifier(yystack_[0].location, yystack_[0].value.as < Symbol > ()); }
#line 982 "tiger_parser.cc"
    break;

  case 23:
#line 167 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new FunCall(yystack_[3].location, yystack_[1].value.as < std::vector<Expr *> > (), yystack_[3].value.as < Symbol > ()); }
#line 988 "tiger_parser.cc"
    break;

  case 24:
#line 171 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new IfThenElse(yystack_[5].location, yystack_[4].value.as < Expr * > (), yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > ()); }
#line 994 "tiger_parser.cc"
    break;

  case 25:
#line 173 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new IfThenElse(yystack_[3].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), new Sequence(nl, {})); }
#line 1000 "tiger_parser.cc"
    break;

  case 26:
#line 178 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, new IntegerLiteral(yystack_[1].location, 0), yystack_[0].value.as < Expr * > (), o_minus); }
#line 1006 "tiger_parser.cc"
    break;

  case 27:
#line 184 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_plus); }
#line 1012 "tiger_parser.cc"
    break;

  case 28:
#line 185 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_minus); }
#line 1018 "tiger_parser.cc"
    break;

  case 29:
#line 186 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_times); }
#line 1024 "tiger_parser.cc"
    break;

  case 30:
#line 187 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_divide); }
#line 1030 "tiger_parser.cc"
    break;

  case 31:
#line 188 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_eq); }
#line 1036 "tiger_parser.cc"
    break;

  case 32:
#line 189 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_neq); }
#line 1042 "tiger_parser.cc"
    break;

  case 33:
#line 190 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_lt); }
#line 1048 "tiger_parser.cc"
    break;

  case 34:
#line 191 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_gt); }
#line 1054 "tiger_parser.cc"
    break;

  case 35:
#line 192 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_le); }
#line 1060 "tiger_parser.cc"
    break;

  case 36:
#line 193 "tiger_parser.yy"
                         { yylhs.value.as < Expr * > () = new BinaryOperator(yystack_[1].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > (), o_ge); }
#line 1066 "tiger_parser.cc"
    break;

  case 37:
#line 194 "tiger_parser.yy"
                         {
        yylhs.value.as < Expr * > () = new IfThenElse(yystack_[1].location, yystack_[2].value.as < Expr * > (),
                            new IfThenElse(yystack_[0].location, yystack_[0].value.as < Expr * > (), new IntegerLiteral(nl, 1), new IntegerLiteral(nl, 0)),
                            new IntegerLiteral(nl, 0));
      }
#line 1076 "tiger_parser.cc"
    break;

  case 38:
#line 199 "tiger_parser.yy"
                        {
        yylhs.value.as < Expr * > () = new IfThenElse(yystack_[1].location, yystack_[2].value.as < Expr * > (),
                            new IntegerLiteral(nl, 1),
                            new IfThenElse(yystack_[0].location, yystack_[0].value.as < Expr * > (), new IntegerLiteral(nl, 1), new IntegerLiteral(nl, 0)));
      }
#line 1086 "tiger_parser.cc"
    break;

  case 39:
#line 208 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new Assign(yystack_[1].location, new Identifier(yystack_[2].location, yystack_[2].value.as < Symbol > ()), yystack_[0].value.as < Expr * > ()); }
#line 1092 "tiger_parser.cc"
    break;

  case 40:
#line 211 "tiger_parser.yy"
                              { yylhs.value.as < Expr * > () = new WhileLoop(yystack_[3].location, yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > ()); }
#line 1098 "tiger_parser.cc"
    break;

  case 41:
#line 215 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new ForLoop(yystack_[7].location, new VarDecl(yystack_[6].location, yystack_[6].value.as < Symbol > (), boost::none, yystack_[4].value.as < Expr * > (), true), yystack_[2].value.as < Expr * > (), yystack_[0].value.as < Expr * > ()); }
#line 1104 "tiger_parser.cc"
    break;

  case 42:
#line 218 "tiger_parser.yy"
                 { yylhs.value.as < Expr * > () = new Break(yystack_[0].location); }
#line 1110 "tiger_parser.cc"
    break;

  case 43:
#line 222 "tiger_parser.yy"
  { yylhs.value.as < Expr * > () = new Let(yystack_[4].location, yystack_[3].value.as < std::vector<Decl *> > (), new Sequence(nl, yystack_[1].value.as < std::vector<Expr *> > ())); }
#line 1116 "tiger_parser.cc"
    break;

  case 44:
#line 225 "tiger_parser.yy"
                              { yylhs.value.as < Expr * > () = new Sequence(yystack_[2].location, yystack_[1].value.as < std::vector<Expr *> > ()); }
#line 1122 "tiger_parser.cc"
    break;

  case 45:
#line 228 "tiger_parser.yy"
       { yylhs.value.as < std::vector<Expr *> > () = std::vector<Expr *>(); }
#line 1128 "tiger_parser.cc"
    break;

  case 46:
#line 229 "tiger_parser.yy"
                  { yylhs.value.as < std::vector<Expr *> > () = yystack_[0].value.as < std::vector<Expr *> > (); }
#line 1134 "tiger_parser.cc"
    break;

  case 47:
#line 232 "tiger_parser.yy"
                    { yylhs.value.as < std::vector<Expr *> > () = std::vector<Expr *>({yystack_[0].value.as < Expr * > ()}); }
#line 1140 "tiger_parser.cc"
    break;

  case 48:
#line 234 "tiger_parser.yy"
  {
    yylhs.value.as < std::vector<Expr *> > () = std::move(yystack_[2].value.as < std::vector<Expr *> > ());
    yylhs.value.as < std::vector<Expr *> > ().push_back(yystack_[0].value.as < Expr * > ());
  }
#line 1149 "tiger_parser.cc"
    break;

  case 49:
#line 240 "tiger_parser.yy"
           { yylhs.value.as < std::vector<Expr *> > () = std::vector<Expr *>(); }
#line 1155 "tiger_parser.cc"
    break;

  case 50:
#line 241 "tiger_parser.yy"
                      { yylhs.value.as < std::vector<Expr *> > () = yystack_[0].value.as < std::vector<Expr *> > (); }
#line 1161 "tiger_parser.cc"
    break;

  case 51:
#line 244 "tiger_parser.yy"
                        { yylhs.value.as < std::vector<Expr *> > () = std::vector<Expr *>({yystack_[0].value.as < Expr * > ()}); }
#line 1167 "tiger_parser.cc"
    break;

  case 52:
#line 246 "tiger_parser.yy"
  {
    yylhs.value.as < std::vector<Expr *> > () = std::move(yystack_[2].value.as < std::vector<Expr *> > ());
    yylhs.value.as < std::vector<Expr *> > ().push_back(yystack_[0].value.as < Expr * > ());
  }
#line 1176 "tiger_parser.cc"
    break;

  case 53:
#line 252 "tiger_parser.yy"
        { yylhs.value.as < std::vector<VarDecl *> > () = std::vector<VarDecl *>(); }
#line 1182 "tiger_parser.cc"
    break;

  case 54:
#line 253 "tiger_parser.yy"
                   { yylhs.value.as < std::vector<VarDecl *> > () = yystack_[0].value.as < std::vector<VarDecl *> > (); }
#line 1188 "tiger_parser.cc"
    break;

  case 55:
#line 256 "tiger_parser.yy"
                      { yylhs.value.as < std::vector<VarDecl *> > () = std::vector<VarDecl *>({yystack_[0].value.as < VarDecl * > ()}); }
#line 1194 "tiger_parser.cc"
    break;

  case 56:
#line 258 "tiger_parser.yy"
  {
    yylhs.value.as < std::vector<VarDecl *> > () = std::move(yystack_[2].value.as < std::vector<VarDecl *> > ());
    yylhs.value.as < std::vector<VarDecl *> > ().push_back(yystack_[0].value.as < VarDecl * > ());
  }
#line 1203 "tiger_parser.cc"
    break;

  case 57:
#line 264 "tiger_parser.yy"
       { yylhs.value.as < std::vector<Decl *> > () = std::vector<Decl *>();}
#line 1209 "tiger_parser.cc"
    break;

  case 58:
#line 266 "tiger_parser.yy"
  {
    yylhs.value.as < std::vector<Decl *> > () = std::move(yystack_[1].value.as < std::vector<Decl *> > ());
    yylhs.value.as < std::vector<Decl *> > ().push_back(yystack_[0].value.as < Decl * > ());
  }
#line 1218 "tiger_parser.cc"
    break;

  case 59:
#line 272 "tiger_parser.yy"
                   { yylhs.value.as < VarDecl * > () = new VarDecl(yystack_[2].location, yystack_[2].value.as < Symbol > (), yystack_[0].value.as < Symbol > (), nullptr); }
#line 1224 "tiger_parser.cc"
    break;

  case 60:
#line 275 "tiger_parser.yy"
                { yylhs.value.as < boost::optional<Symbol> > () = boost::none; }
#line 1230 "tiger_parser.cc"
    break;

  case 61:
#line 276 "tiger_parser.yy"
             { yylhs.value.as < boost::optional<Symbol> > () = yystack_[0].value.as < Symbol > (); }
#line 1236 "tiger_parser.cc"
    break;


#line 1240 "tiger_parser.cc"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[+yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yy_error_token_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yy_error_token_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
    goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
  `-------------------------------------*/
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;


  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
  `-----------------------------------*/
  yyabortlab:
    yyresult = 1;
    goto yyreturn;


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  tiger_parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

  // Generate an error message.
  std::string
  tiger_parser::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    std::ptrdiff_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state merging
         (from LALR or IELR) and default reductions corrupt the expected
         token list.  However, the list is correct for canonical LR with
         one exception: it will still contain any token that will not be
         accepted due to an error action in a later state.
    */
    if (!yyla.empty ())
      {
        symbol_number_type yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];

        int yyn = yypact_[+yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yy_error_token_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
      default: // Avoid compiler warnings.
        YYCASE_ (0, YY_("syntax error"));
        YYCASE_ (1, YY_("syntax error, unexpected %s"));
        YYCASE_ (2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_ (3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_ (4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_ (5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    std::ptrdiff_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char tiger_parser::yypact_ninf_ = -37;

  const signed char tiger_parser::yytable_ninf_ = -1;

  const short
  tiger_parser::yypact_[] =
  {
      43,    43,    43,    43,    43,   -36,   -37,   -37,    11,   -37,
     -37,    23,   150,   -37,   -37,   -37,   -37,   -37,   -37,   -37,
     -37,   -37,   -37,   -37,   -37,   -37,   150,    18,    21,   -37,
     135,    -5,     5,   -13,    43,    43,   -37,    43,    43,    43,
      43,    43,    43,    43,    43,    43,    43,    43,    43,   -37,
      43,    43,    43,    43,    43,    -9,    -8,   -37,   -37,   -37,
     150,    24,    27,   150,     7,     7,   -37,   -37,   162,   162,
     162,   162,   162,   162,   172,    45,   150,   119,   150,   100,
       4,    26,    63,   -37,    43,    43,    43,   -37,    31,    34,
      50,   150,   150,    77,    70,    68,    74,   -37,   -37,    43,
      43,    41,    63,    31,   150,   150,   -37,    65,   -37,    43,
     150
  };

  const signed char
  tiger_parser::yydefact_[] =
  {
       0,    45,     0,     0,     0,     0,    57,    42,    22,    20,
      21,     0,     2,     5,     6,     8,     9,    12,    11,    10,
      13,    14,    15,    16,    17,     7,    47,     0,    46,    26,
       0,     0,     0,     0,    49,     0,     1,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    44,
       0,     0,     0,     0,    45,     0,     0,    58,     3,     4,
      51,     0,    50,    39,    27,    28,    29,    30,    31,    32,
      33,    35,    34,    36,    37,    38,    48,    25,    40,     0,
       0,     0,    60,    23,     0,     0,     0,    43,    53,     0,
       0,    52,    24,     0,     0,     0,    54,    55,    61,     0,
       0,     0,    60,     0,    18,    41,    59,     0,    56,     0,
      19
  };

  const signed char
  tiger_parser::yypgoto_[] =
  {
     -37,   -37,   -37,     0,   -37,   -37,   -37,   -37,   -37,   -37,
     -37,   -37,   -37,   -37,   -37,   -37,   -37,   -37,   -37,    29,
     -37,   -37,   -37,   -37,   -37,   -37,    -2,     1
  };

  const signed char
  tiger_parser::yydefgoto_[] =
  {
      -1,    11,    57,    26,    58,    59,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    27,
      28,    61,    62,    95,    96,    33,    97,    90
  };

  const signed char
  tiger_parser::yytable_[] =
  {
      12,    32,    29,    30,    31,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    34,    54,    39,
      40,    55,    56,    36,    52,    49,    50,    53,    81,    82,
      84,    83,    88,    35,    60,    63,    87,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    75,     1,
      76,    77,    78,    79,     2,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,     3,    89,    94,     4,
       5,    98,    99,     6,   101,   102,     7,   103,   106,   109,
       8,     9,    10,    80,    91,    92,    93,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,   104,
     105,   108,     0,   107,     0,     0,   100,     0,     0,   110,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,     0,     0,     0,     0,     0,     0,    86,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,     0,     0,     0,    85,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,     0,     0,    51,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    37,    38,    39,    40,    -1,    -1,    -1,    -1,
      -1,    -1,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46
  };

  const signed char
  tiger_parser::yycheck_[] =
  {
       0,    37,     2,     3,     4,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,     6,    31,    12,
      13,    34,    35,     0,    29,     7,     5,    22,    37,    37,
       3,     7,     6,    22,    34,    35,    32,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,     6,
      50,    51,    52,    53,    11,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    23,     4,    37,    26,
      27,    37,    22,    30,     4,     7,    33,     3,    37,    14,
      37,    38,    39,    54,    84,    85,    86,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    99,
     100,   103,    -1,   102,    -1,    -1,    29,    -1,    -1,   109,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    -1,    -1,    -1,    -1,    -1,    -1,    28,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    -1,    -1,    -1,    25,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    -1,    -1,    24,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19
  };

  const signed char
  tiger_parser::yystos_[] =
  {
       0,     6,    11,    23,    26,    27,    30,    33,    37,    38,
      39,    43,    45,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    45,    61,    62,    45,
      45,    45,    37,    67,     6,    22,     0,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,     7,
       5,    24,    29,    22,    31,    34,    35,    44,    46,    47,
      45,    63,    64,    45,    45,    45,    45,    45,    45,    45,
      45,    45,    45,    45,    45,    45,    45,    45,    45,    45,
      61,    37,    37,     7,     3,    25,    28,    32,     6,     4,
      69,    45,    45,    45,    37,    65,    66,    68,    37,    22,
      29,     4,     7,     3,    45,    45,    37,    69,    68,    14,
      45
  };

  const signed char
  tiger_parser::yyr1_[] =
  {
       0,    42,    43,    44,    44,    45,    45,    45,    45,    45,
      45,    45,    45,    45,    45,    45,    45,    45,    46,    47,
      48,    49,    50,    51,    52,    52,    53,    54,    54,    54,
      54,    54,    54,    54,    54,    54,    54,    54,    54,    55,
      56,    57,    58,    59,    60,    61,    61,    62,    62,    63,
      63,    64,    64,    65,    65,    66,    66,    67,    67,    68,
      69,    69
  };

  const signed char
  tiger_parser::yyr2_[] =
  {
       0,     2,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     5,     8,
       1,     1,     1,     4,     6,     4,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       4,     8,     1,     5,     3,     0,     1,     1,     3,     0,
       1,     1,     3,     0,     1,     1,     3,     0,     2,     3,
       0,     2
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const tiger_parser::yytname_[] =
  {
  "\"end of file\"", "error", "$undefined", "\",\"", "\":\"", "\";\"",
  "\"(\"", "\")\"", "\"{\"", "\"}\"", "\"+\"", "\"-\"", "\"*\"", "\"/\"",
  "\"=\"", "\"<>\"", "\"<\"", "\"<=\"", "\">\"", "\">=\"", "\"&\"",
  "\"|\"", "\":=\"", "\"if\"", "\"then\"", "\"else\"", "\"while\"",
  "\"for\"", "\"to\"", "\"do\"", "\"let\"", "\"in\"", "\"end\"",
  "\"break\"", "\"function\"", "\"var\"", "\"uminus\"", "\"id\"",
  "\"string\"", "\"int\"", "TYPE", "OF", "$accept", "program", "decl",
  "expr", "varDecl", "funcDecl", "stringExpr", "intExpr", "var",
  "callExpr", "ifExpr", "negExpr", "opExpr", "assignExpr", "whileExpr",
  "forExpr", "breakExpr", "letExpr", "seqExpr", "exprs", "nonemptyexprs",
  "arguments", "nonemptyarguments", "params", "nonemptyparams", "decls",
  "param", "typeannotation", YY_NULLPTR
  };

#if YYDEBUG
  const short
  tiger_parser::yyrline_[] =
  {
       0,   122,   122,   125,   126,   129,   130,   131,   132,   133,
     134,   135,   136,   137,   138,   139,   140,   141,   144,   148,
     154,   158,   162,   166,   170,   172,   177,   184,   185,   186,
     187,   188,   189,   190,   191,   192,   193,   194,   199,   207,
     211,   214,   218,   221,   225,   228,   229,   232,   233,   240,
     241,   244,   245,   252,   253,   256,   257,   264,   265,   272,
     275,   276
  };

  // Print the state stack on the debug stream.
  void
  tiger_parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  tiger_parser::yy_reduce_print_ (int yyrule)
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG


} // yy
#line 1721 "tiger_parser.cc"

#line 279 "tiger_parser.yy"


void
yy::tiger_parser::error (const location_type& l,
                          const std::string& m)
{
  utils::error (l, m);
}
