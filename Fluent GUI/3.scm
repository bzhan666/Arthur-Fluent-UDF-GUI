(define apply-cb #t)
(define update-cb #f)

(define table)
(define myList1)
(define myList2)

(define (button-cb . args)
	(cx-set-list-items myList2(cx-show-list-selections myList1))
	)
(define bzhan666 (cx-create-panel "BZHAN666" apply-cb update-cb))
(set! table (cx-create-table bzhan666 "CFD"))
(set! myList1 (cx-create-drop-down-list table "List 1"  'multiple-selections #t 'row 0))
(cx-set-list-items myList1 (list "EDC1" "EDC2" "EDC3" "EDC4" "EDC5"))

(set! myList2 (cx-create-drop-down-list table "List 2"  'multiple-selections #t 'row 1))

(cx-create-button table "Button" 'activate-callback button-cb 'row 2)
(cx-show-panel bzhan666)
