(define apply-cb #t)
(define update-cb #f)

(define table)
;(define txtField)

;(define counter 0)
(define (button-cb . args)
	(ti-info "nishigedashuaib")
 ;(set! counter (+ counter 1))
 ;(cx-set-text-entry txtField (string-append "Pengyuyan: " (number->string counter)))

)

(define my-dialog-box (cx-create-panel "Bzhan666" apply-cb update-cb))
(set! table (cx-create-table my-dialog-box "CFD"))

;(set! txtField (cx-create-text-entry table "" 'row 0 'col 0))
(cx-create-button table "Button" 'activate-callback button-cb 'row 1 'col 0)

(cx-show-panel my-dialog-box)