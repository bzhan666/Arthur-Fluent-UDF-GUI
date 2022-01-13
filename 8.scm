(define classic
	(lambda ()
		(ti-menu-load-string(format #f "display/set/colors/color-scheme/classic \n"))
		(ti-menu-load-string(format #f "display/set/colors/background \"black\" "))
		(ti-menu-load-string(format #f "display/set/colors/foreground \"white\" "))
		)
	)

(define workbench
	(lambda ()
		(ti-menu-load-string(format #f "display/set/colors/color-scheme/workbench \n"))
		(ti-menu-load-string(format #f "display/set/colors/foreground \"black\" "))
		)
	)

(define whiteColor
	(lambda()
		(ti-menu-load-string(format #f "display/set/colors/color-scheme/classic \n"))
		(ti-menu-load-string(format #f "display/set/colors/background \"white\" "))
		(ti-menu-load-string(format #f "display/set/colors/foreground \"black\" "))
		)
	)

(cx-add-menu "Background" #\y)
(cx-add-item "Background" "Workbench" #f #f #t workbench)
(cx-add-item "Background" "Classic" #f #f #t classic)
(cx-add-item "Background" "White" #f #f #t whiteColor)