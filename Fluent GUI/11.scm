(define (apply-cb) #t)
(define update-cb #f)

(define table)
(define int1)
(define int2)
(define int3)
(define int4)

(define my-dialog-box (cx-create-panel "My Dialog Box" apply-cb update-cb))
(set! table (cx-create-table my-dialog-box "This is an example Dialog Box"))
(set! int1 (cx-create-integer-entry table "Row 1 Col 1" 'row 1 'col 1))
(set! int2 (cx-create-integer-entry table "Row 1 Col 2" 'row 1 'col 2))
(set! int3 (cx-create-integer-entry table "Row 2 Col 1" 'row 2 'col 1))
(set! int4 (cx-create-integer-entry table "No Row/Col Specified"))

(cx-show-panel my-dialog-box)
