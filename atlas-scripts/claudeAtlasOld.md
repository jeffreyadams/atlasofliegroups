# claudeAtlas.md — Generic guide to running the Atlas software

A reusable reference (not tied to any one project) for driving the **atlas**
program — the interpreter for the *axis* language plus the Atlas of Lie Groups
and Representations library. Written for Claude Code sessions, but useful to
anyone.

---

## 1. What atlas is

`atlas` is an expression-oriented interpreter. You type expressions; it type-checks
them, evaluates them, and prints a `Value:` line. On top of a general-purpose
programming language called **axis**, it adds many representation-theory value
types (`RealForm`, `Param`, `RootDatum`, `KGBElt`, `ParamPol`, …) and library
functions operating on them.

The `.at` files in `atlas-scripts/` are axis source files — a large standard
library of user-defined functions layered on top of the built-ins.

---

## 2. Running it

The binary lives at `~/atlasSoftware/master/atlas` and scripts at
`~/atlasSoftware/master/atlas-scripts/`. The standard idiom is to run **from
inside `atlas-scripts/`** so relative `<file.at` includes resolve:

```bash
cd ~/atlasSoftware/master/atlas-scripts
../atlas all.at
```

`all.at` is a manifest that `<`-includes the full recommended library
(`basic.at`, `hermitian.at`, `W_reps.at`, `associated_variety_annihilator.at`,
…). Loading it gives you essentially everything. It takes a little time to load.

**Interactive:** you land at an `atlas>` prompt. Type expressions, end with `quit`.

**Batch / non-interactive (best for scripted or Claude-driven use):** feed
commands on stdin and let it exit with `quit`:

```bash
cd ~/atlasSoftware/master/atlas-scripts
printf 'set G=Sp(4,R)\nset p=trivial(G)\np\ndimension(G)\nquit\n' | ../atlas all.at
```

Typical output:
```
Variable G: RealForm (overriding previous instance, which had type RealForm)
Variable p: Param
Value: final parameter(x=10,lambda=[2,1]/1,nu=[2,1]/1)
Value: 10
Bye.
```

Tips for batch use:
- Always end with `quit` so the process exits.
- Wrap in a `timeout` (e.g. `timeout 300 ../atlas all.at`) — some computations
  (large blocks, E8, cells) can run a very long time.
- To run a whole script file, either pass it after `all.at`
  (`../atlas all.at myscript.at`) or, inside a session, use `<myscript.at`.

**Installed command:** if `make install` was run, a wrapper called `atlas` may
be on `$PATH`; it preloads a standard prelude and can be called from any
directory (`atlas` instead of `../atlas`). When in doubt, use the explicit
`cd atlas-scripts && ../atlas all.at` form above.

---

## 3. Loading files

- `<filename` — read `filename` as a series of commands (the `.at` suffix is
  added automatically if needed). Files may include other files.
- **Include-once:** after a file loads successfully, atlas remembers it and
  silently skips later `<` of the same file — so scripts can freely `<`-include
  their dependencies without re-running them.
- `<<filename` — **force** re-read even if already loaded. Use this after you
  have **edited** a `.at` file and want the new version in the current session.
- Output redirection: prefix an expression (not a `set`) with `>file` (create)
  or `>>file` (append) to send its output to a file:
  `>out.txt block_sizes(ic)`.

---

## 4. The axis language — essentials

### Types
Strongly and **statically typed**: expressions are fully type-checked *before*
any evaluation, so type errors are caught up front. Primitive types include
`bool`, `int` (unbounded), `rat` (exact rationals), `string`, `vec`, `mat`.
Composite types include rows `[T]`, tuples `(T1,T2,…)`, rational vectors
`ratvec`, plus all the Atlas-specific types.

Type errors look like:
```
Type error: ... has wrong type: found RootDatum while InnerClass was needed.
```

### Defining names
- `set name = expr` — define/replace a **global** identifier.
- `let x = e1 in  <expr using x>` — a **local** binding, scoped to the `in` body.
  Chain locals with `then` (= "in let"): `let a=… then b=f(a) in g(a,b)`.
- `:=` — **assignment** to an existing variable (global or local). `=` alone is
  equality; you need `set`/`let` to introduce, `:=` to mutate.
  Compound forms: `x +:= y` means `x := x + y` (works for any operator).
- `!name` — a **constant** binding that forbids later assignment
  (`set !s = Split:(0,1)`).

Value semantics: names never share storage. Assigning to one never affects
another; functions get arguments strictly by value and cannot mutate the
caller's variables — pass results back via return values.

### Functions
```
set f (LieType lt) = void:
   let rd = simply_connected(lt)
   then ic = inner_class(rd,"s")
   in print_block(block(ic, ...))
```
Syntax: `set name (Type arg, Type arg, …) = body`. The optional `void:` (or any
`Type:`) prefix is a **cast** documenting the return type; usually axis infers
it. atlas replies `Defined f: (LieType->)` telling you the inferred signature.

### Overloading
Almost every function/operator is **overloaded** — the same name can have many
type signatures, and your `set` adds a signature rather than clobbering the
built-in. Implicit conversions (e.g. `string`→`LieType`, `int`→`rat`) are
inserted automatically where needed.

### Control structures
- Sequencing with `;` (value is the last expr); `next` returns the *first*
  operand's value. `begin … end` groups like parentheses.
- `if cond then … elif cond then … else … fi`
- Loops: `while cond do … od`, `for x in list do … od`, `for i:n do … od`.
  Unusually, loops **collect** their body values into a row `[…]`.

---

## 5. Handy REPL commands

| Command | Effect |
|---|---|
| `quit` | exit (`Bye.`) |
| `whattype expr` | type-check and print the type — no evaluation |
| `whattype name ?` | list all overloads of a function/operator name |
| `showall` | dump type + value of every defined identifier |
| `set verbose` / `set quiet` | toggle verbose type-checker output |
| `>file expr` / `>>file expr` | redirect an expression's output to a file |
| `<file` / `<<file` | read commands from a file / force re-read |

Example:
```
atlas> whattype + ?
Overloaded instances of +
  (int,int)->int
  (ratvec,ratvec)->ratvec
  (ParamPol,Param)->ParamPol
  ...
```

---

## 6. Conventions & gotchas

- Functions named `print_…` / `show_…` **print** and return `void` — their
  output is not a `Value:` you can capture or reuse. To manipulate a result,
  use the non-printing function that returns the actual value.
- A `Param` prints as `final parameter(x=…, lambda=…, nu=…)`; "final" means it
  is in the canonical form the library uses for irreducibles.
- Because type-checking precedes evaluation, a whole batch fails at parse/type
  time if any command is ill-typed — check early commands first when debugging.
- An error while reading an included file **abandons all active input files**
  and returns you to the prompt; the partially-read file is *not* marked loaded,
  so a plain `<` will re-read it.
- Building groups: `Sp(4,R)`, `SL(2,R)`, `GL(n)`, `inner_class("B5.A3.T1", …)`,
  `simply_connected(lt)`, etc. `rho(G)`, `trivial(G)`, `dimension(G)` are common
  starting points.

---

## 7. Reference material on disk

Inside `atlas-scripts/`:
- `atlas.help` — the axis language reference + atlas intro (the source for this
  guide; sections: Types, Identifiers, Functions, Control structures, I/O).
- `atlas-functions.help` — catalog of library functions.
- `new-types.help` — the Atlas-specific value types.
- `messages/*.help` — per-topic help (block, KGB, KL, branch, …).

To find where a function is defined or how it's used:
```bash
grep -rn "function_name" ~/atlasSoftware/master/atlas-scripts/*.at
```
