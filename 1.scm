(define apply-cb #t)
(define update-cb #f)

(define bzhan666 (cx-create-panel "bzhan666" apply-cb update-cb))
(define table)
(define checkbox)
(define ints)
(define reals)
(define txt)

(set! table(cx-create-table bzhan666 "cfd"))

(set! checkbox(cx-create-toggle-button table "Pengyuyan" 'row 0))
(set! ints(cx-create-integer-entry table "Wuyanzhu" 'row 1))
(set! reals(cx-create-real-entry table "EDC" 'row 2))
(set! txt (cx-create-text-entry table "For you" 'row 3))

(cx-set-integer-entry ints 1)
(cx-set-real-entry reals 1.2)
(cx-set-text-entry txt "da shuai b")
(cx-show-panel bzhan666)