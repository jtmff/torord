<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.2/jquery.min.js"></script>
<script type="text/javascript">
 function changeImage(current) {
	var imagesNumber = 20;

	for (i=1; i<=imagesNumber; i++) {
		if (i == current) {
			document.getElementById("normal" + current).style.display = "block";
		} else {
			document.getElementById("normal" + i).style.display = "none";
		}
	}
}



$(document).ready(function(){
    /*$("p").click(function(){
        $(this).hide();
    });*/
	$('input').on('change', function () {
		var v = $(this).val();
		//$(".jqTest").css('font-size', v + 'px')
		//$('span').html(v);
		$('#thumbs img').width(10*v); // Units are assumed to be pixels
		//$('#thumbs').height(700);
	});
	
	$('#imageContainer img').each(function (index) {
		if ($(this).attr('onclick') != null) {                    
			if ($(this).attr('onclick').indexOf("runThis()") == -1) {                        
				$(this).click(function () {
					$(this).attr('onclick');
					var src = $(this).attr("src");
					var bigImage = $(this).attr("openedImage");
					//$('span').html(bigImage);
					if (typeof bigImage !== typeof undefined && bigImage !== false) {
						if (bigImage == "src_openedImage") {
							ShowLargeImage2(src);
						} else {
							ShowLargeImage2(bigImage);
						}
					} /*else {
						ShowLargeImage(src);					
					}*/
				});
			}
		} else {                    
			$(this).click(function () {                        
				var src = $(this).attr("src");
				var bigImage = $(this).attr("openedImage");
				//$('span').html(bigImage);
				if (typeof bigImage !== typeof undefined && bigImage !== false) {
					if (bigImage == "src_openedImage") {
						ShowLargeImage2(src);
					} else {
						ShowLargeImage2(bigImage);
					}
				} /*else {
					ShowLargeImage(src);					
				}*/
			});
		}
	});

	$('body').on('click', '.modal-overlay', function () {
		$('.modal-overlay, .modal-img').remove();
	});

	function ShowLargeImage(imagePath) {
		$('body').append('<div class="modal-overlay"></div><div class="modal-img"><img src="' + imagePath.replace("small","large") + '" /></div>'); 
		$('.modal-img').animate({
			opacity: 1
		},1000);
		$('.modal-img').click(function () {
			$('.modal-overlay, .modal-img').remove();
		});
	}
	function ShowLargeImage2(imagePathBig) {
		$('body').append('<div class="modal-overlay"></div><div class="modal-img"><img src="' + imagePathBig + '" /></div>'); 
		$('.modal-img').animate({
			opacity: 1
		},1000);
		$('.modal-img').click(function () {
			$('.modal-overlay, .modal-img').remove();
		});
	}
});
</script>