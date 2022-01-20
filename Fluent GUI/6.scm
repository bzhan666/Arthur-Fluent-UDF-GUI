(define apply-cb #t)
(define update-cb #f)

(define table)
(define myDropList1)
(define myDropList2)

(define (button-cb . args)
	(cx-set-list-items myDropList2 (cx-show-list-selections myDropList1))
)

(define bzhan666 (cx-create-panel "BZHAN666" apply-cb update-cb))
(set! table (cx-create-table bzhan666 "CFD"))

(set! myDropList1 (cx-create-drop-down-list table "Drop Down List 1" 'row 0))
(cx-set-list-items myDropList1 (list "EDC1" "EDC2" "EDC3" "EDC4" "EDC5"))

(set! myDropList2 (cx-create-drop-down-list table "Drop Down List 2" 'row 1))

(cx-create-button table "Button" 'activate-callback button-cb 'row 2)

(cx-show-panel bzhan666)
