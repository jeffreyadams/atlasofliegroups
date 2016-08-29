;;; atlas.el  --- Support for editing atlas scripts

;; Author Marc van Leeuwen <Marc.van-Leeuwen@math.univ-poitiers.fr>


(defconst atlas-font-lock-externals
 (let ()
  (list
   ;; include file directives
   (list (rx
	  line-start (zero-or-more blank) (repeat 1 2 "<") (zero-or-more blank)
	  (group (one-or-more (not whitespace)))
	  ) 1 'font-lock-preprocessor-face)
   ;'("^\\s-*<<?\\(\\S-+\\)" 1 font-lock-preprocessor-face)

   ;; type definitions
   (list (rx
	  line-start (zero-or-more blank) ":" (zero-or-more blank)
	  (group (one-or-more (or wordchar (syntax symbol))))
	  ) 1 'font-lock-type-face)
   ;'("^\\s-*:\\s-+\\(\\(\\w\\|\\s_\\)+\\)" 1 font-lock-preprocessor-face)

   ; function definitions
   (list (rx
	  line-start (zero-or-more blank) "set" (one-or-more blank)
	  (group (one-or-more (or wordchar (syntax symbol))))
	  ) 1 'font-lock-function-name-face)
   ; '("^\\s-*set\\s-+\\(\\(\\w\\|\\s_\\)+\\)" 1 font-lock-function-name-face)
   )
  )
 "Matching structure for axis outer level structures"
)

(defconst atlas-font-lock-keywords
  (let ((keywords
	 '("quit" "set" "let" "in" "begin" "end"
	   "if" "then" "else" "elif" "fi" "and" "or" "not"
	   "next" "do" "dont" "from" "downto" "while" "for" "od"
	   "case" "esac" "rec_fun" "true" "false"  "die"
	   "break" "return" "whattype" "showall" "forget"))
	(types '("int" "rat" "string" "bool" "void"
		 "vec" "mat" "ratvec"
		 "LieType" "RootDatum" "InnerClass" "RealForm"
		 "CartanClass" "KGBElt" "Block" "Param" "Split" "ParamPol"))
	)
    (list (regexp-opt keywords 'words)
	  (cons (regexp-opt types 'words) 'font-lock-type-face)
	  (cons "~" 'font-lock-negation-char-face)
	  )
   )
  "Matching structure for altas keywords and predefined types"
)

(define-derived-mode atlas-mode fundamental-mode "Atlas"
 "Major mode for editing Atlas scripts"
 (let ((syntax-alist
	'( ( ?\{ . "<}n" ) (?\} . ">{n")
	   (?\" . "\"")
	   (?_ . "_")
           (?= . ".")
	   ))
       (level1 (append atlas-font-lock-externals atlas-font-lock-keywords))
       (level2 (append atlas-font-lock-externals atlas-font-lock-keywords))
       )
  (setq font-lock-defaults (list
			    (list 'atlas-font-lock-externals level1 level2)
			    nil nil syntax-alist)))
)